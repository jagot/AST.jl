using AtomicLevels

function nonh()
    info("Executing nonh")
    # Angular integrals
    clean_mpi_tmp(r"fort.[0-9]+",
                  r"yint.lst.[0-9]+",
                  r"c.lst.[0-9]+")
    mpi_run(mpi_cmd("$atsp/nonh_mpi"))
    info("Finished executing nonh")
end

function hf(name, term::Term, Z,
            closed::Config, ref_set::Config,
            vary = "all",
            default_electron_parameters = true,
            rem_parameters = ["y","n","T,1.D-8","y","y"],
            additional_parameters = [])
    closed = join([@sprintf("%2d%s", c[1], ells[c[2]+1]) for c in closed], " ")

    pipe_file_run("$atsp/hf",
                  """$name,$(string(term,false)),$(float(Z)) ! Name, final term, Z (next row: closed orbitals)
 $closed
$(orbital_string(ref_set, true, false; escape_principal_qn = true)) ! Electrons outside closed orbitals
$vary ! Which orbitals to vary
$(y_or_n(default_electron_parameters)) ! Default electron parameters
$(y_or_n(length(rem_parameters)==0)) ! Remaining parameters
$(join(rem_parameters, "\n"))
$(y_or_n(length(additional_parameters)>0)) ! Additional parameters
$(join(additional_parameters, "\n"))
n ! Don't continue along the sequence
""")
end

function mchf(name,Z,
              term_weights,
              active::Config,
              spectroscopic::Config = [],
              iter = 200,
              cfg_tol = 1e-8,
              scf_tol = 1e-8)
    pipe_file_run("$atsp/mchf_mpi",
                  """$name,$(float(Z))
$(join(term_weights, '\n'))
$(orbital_string(active,false,false,","; escape_principal_qn = true))
$(orbital_string(spectroscopic,false,false,"," ; escape_principal_qn = true))
y ! Default electron parameters
n ! Custom other parameters
y ! Default values for NO,REL,STRONG
n ! Custom values for PRINT,CFGTOL,SCFTOL
F,$cfg_tol,$scf_tol
n ! Custom values for NSCF,IC
$iter,0
y ! Default values for ACFG,LD,TRACE
""")
end

function atsp_save_run(config, term)
    cpf("wfn.out", "wfn.inp")
    name = active_file(config, term)
    cpf("wfn.out", "$name.w")
    cpf("summry", "$name.s")
    cpf("clist.out", "$name.c")
    cpf(filter(f -> ismatch(Regex("^$(string(term)).*\\.l\$"), f), readdir())[1], "$name.l")
end

function breit_pauli(config, term;
                     restarting = false,
                     guesses = false,
                     which_op = 2,
                     all_rel = true,
                     all_inter = true,
                     atom_mass = Inf)
    name = active_file(config, term)

    (isfile("$name.l") && isfile("$name.c")) || error("Could not find $name.l/c")
    cpf("$name.l", "$name.l.bak")
    Js = J_range(term)
    j2 = map(Js) do JJ
        round(Int, 2JJ)
    end
    eiv = join([1 for JJ in reverse(Js)], '\n')
    println("$term ⟹ J ∈ $Js")

    atom_mass_s = (atom_mass == Inf) ? "" : "\n$atom_mass"

    pipe_file_run("$atsp/bpci",
                  """$name, y, n ! Atom, relativistic, mass correction
$(y_or_n(restarting)) ! Restarting
$(y_or_n(guesses)) ! Existing $name.l/.j as guesses
$(maximum(j2)),$(minimum(j2)) ! max(2J),min(2J)
$eiv ! Eigenvalues
$which_op ! Non-relativistic operators and selected relativistic
$(y_or_n(all_rel)) ! All relativistic operators
$(y_or_n(all_inter)) ! All interactions
$(y_or_n(atom_mass == Inf)) ! Default Rydberg constant$(atom_mass_s)
""")
    cpf("$name.l.bak", "$name.l")
end

function atsp_cp_wfn(a,b)
    mkpath("$b")
    cpf("$a/wfn.out", "$b/wfn.inp")
end

function read_hf_eng()
    # Read the total energy from hf.log
    open("hf.log") do file
        float(split(readall(file))[end-11])
    end
end

function read_eng_blocks(f::Function, config, term, ending)
    open("$(active_file(config, term)).$ending") do file
        f() do
            while !ismatch(r"^[ ]*2\*J =", readline(file)); end
            readline(file)
            float(split(readline(file))[2])
        end
    end
end

function read_mchf_eng(config, term)
    # Read the total (non-relativistic) energy from $name.l
    read_eng_blocks(config, term, "l") do read_block
        read_block()
    end
end

function read_breit_pauli_eng(config, term)
    Js = J_range(term)
    # Read the total (Breit–Pauli corrected) energies from $name.j
    read_eng_blocks(config, term, "j") do read_block
        [read_block()
         for JJ in reverse(Js)]
    end
end

function atsp_clean()
    # Clean up temporary files that take up a lot of space
    for f in filter(f -> contains(f, ".lst"), readdir())
        rm(f)
    end
end


function hf_mchf_bp(config::Config,
                    term::Term,
                    Z,
                    ncorr,
                    wfn0 = nothing;
                    active::Function = active_set,
                    overwrite = true,
                    atom_mass = Inf,
                    csf_filter = (a...) -> true,
                    hf_core = true,
                    bp = true)
    println(repeat("=", 80))
    conf = active_file(config, term)
    println(config, " ", term)

    overwrite_i = i -> overwrite && isfile("$i/wfn.out") || !isfile("$i/wfn.out")
    dir_run(conf) do
        energies = zeros(ncorr+1)

        Js = J_range(term)
        bp_energies = zeros(ncorr+1, length(Js))

        # Hartree--Fock run
        overwrite_i(0) && dir_run("0") do
            if wfn0 != nothing
                println("Initial guess: $wfn0")
                cpf(wfn0, "wfn.inp")
            end
            csfgenerate(active(config, 0),
                        csf_filter,
                        (term, 0, config))
            nonh()
            hf(conf, term, Z,
               closed(config),
               open(config),
               hf_core ? "all" : "=$(length(open(config)))")
            atsp_clean()
        end

        dir_run("0") do
            energies[1] = read_hf_eng()
            bp_energies[1,:] = energies[1]
        end

        # Multiconfigurational Hartree--Fock, layer by layer, with
        # additional Breit--Pauli calculation at the end of each run
        for i = 1:ncorr
            atsp_cp_wfn(i-1,i)
            overwrite_i(i) && dir_run("$i") do
                act = csfgenerate(active(config, i),
                                  csf_filter,
                                  (term, 2, config))
                nonh()
                mchf(conf, Z, 1,
                     act, config)
                atsp_save_run(config, term)
                bp && breit_pauli(config, term; atom_mass = atom_mass)
                atsp_clean()
            end

            dir_run("$i") do
                energies[i+1] = read_mchf_eng(config, term)
                bp && (bp_energies[i+1,:] = read_breit_pauli_eng(config, term))
            end
        end

        save_eng([energies bp_energies]')
    end
    conf, abspath("$conf/$ncorr/wfn.out")
end

export nonh, hf, mchf,
atsp_save_run, breit_pauli,
atsp_cp_wfn, read_hf_eng, read_mchf_eng, read_breit_pauli_eng,
atsp_clean,
hf_mchf_bp
