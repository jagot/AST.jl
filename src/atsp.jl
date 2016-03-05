function csfgenerate_input(active,
                           configurations,
                           term,
                           nexc)
    """$(join(configurations, '\n')) ! Configurations

$active
$term ! Resulting term
$nexc ! Number of excitations"""
end

function csfgenerate(active, lists)
    # If core not in list, the value is zero (no core)
    lists = join(map(l -> csfgenerate_input(active, l...), lists), "\ny\n")
    pipe_file_run("$atsp/csfexcitation",
                  """$lists
n ! No more lists

""")
    cpf("csfexcitation.log", "$(active_file(active)).exc")
    lsgen_inp = split(readall(open("excitationdata.sh")), "\n")
    pipe_file_run("$atsp/lsgen",
                  join(lsgen_inp[3:end-2], "\n"))
    cpf("clist.out", "cfg.inp")
end

function nonh()
    # Angular integrals
    clean_mpi_tmp(r"fort.[0-9]+",
                  r"yint.lst.[0-9]+",
                  r"c.lst.[0-9]+")
    mpi_run(mpi_cmd("$atsp/nonh_mpi"))
end

function hf(name, term, Z,
            closed, ref_set,
            vary = "all",
            default_electron_parameters = true,
            default_rem_parameters = true,
            additional_parameters = [])
    closed = join([@sprintf("%3s", c) for c in closed], " ")

    pipe_file_run("$atsp/hf",
                  """$name,$term,$(float(Z)) ! Name, final term, Z (next row: closed orbitals)
 $closed
$ref_set ! Electrons outside closed orbitals
$vary ! Vary $all orbitals
$(y_or_n(default_electron_parameters)) ! Default electron parameters
$(y_or_n(default_rem_parameters)) ! Default values for remaining parameters
$(y_or_n(length(additional_parameters)>0)) ! Additional parameters
$(join(additional_parameters, "\n"))
n ! Don't continue along the sequence
""")
end

function mchf(name,Z,
              term_weights,
              active,
              spectroscopic = "",
              iter = 200,
              cfg_tol = 1e-8,
              scf_tol = 1e-8)
    pipe_file_run("$atsp/mchf_mpi",
                  """$name,$(float(Z))
$(join(term_weights, '\n'))
$active
$spectroscopic
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

function atsp_save_run(active, term)
    cpf("wfn.out", "wfn.inp")
    active = active_file(active)
    cpf("wfn.out", "$active.w")
    cpf("summry", "$active.s")
    cpf("clist.out", "$active.c")
    cpf(filter(f -> ismatch(Regex("^$term.\\.l\$"), f), readdir())[1], "$active.l")
end

function breit_pauli(active, term,
                     restarting = false,
                     guesses = false,
                     which_op = 2,
                     all_rel = true,
                     all_inter = true,
                     def_Rydberg = true)
    active = active_file(active)

    println(pwd(), " ", isfile("$active.l"))
    isfile("$active.l") || error("Could not find $active.l")
    j2 = term_to_2j_range(term)
    eiv = join([1 for jj in maximum(j2):-2:minimum(j2)], '\n')
    println("$term => J ∈ $(j2_to_jstr(j2))")

    pipe_file_run("$atsp/bpci",
                  """$active, y, n ! Atom, relativistic, mass correction
$(y_or_n(restarting)) ! Restarting
$(y_or_n(guesses)) ! Existing $active.l/.j as guesses
$(maximum(j2)),$(minimum(j2)) ! max(2J),min(2J)
$eiv ! Eigenvalues
$which_op ! Non-relativistic operators and selected relativistic
$(y_or_n(all_rel)) ! All relativistic operators
$(y_or_n(all_inter)) ! All interactions
$(y_or_n(def_Rydberg)) ! Default Rydberg constant
""")
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

function read_eng_blocks(f::Function, active, ending)
    open("$(active_file(active)).$ending") do file
        f() do
            while !ismatch(r"^[ ]*2\*J =", readline(file)); end
            readline(file)
            float(split(readline(file))[2])
        end
    end
end

function read_mchf_eng(active)
    # Read the total (non-relativistic) energy from $active.l
    read_eng_blocks(active, "l") do read_block
        read_block()
    end
end

function read_breit_pauli_eng(active, term)
    j2 = term_to_2j_range(term)
    # Read the total (Breit–Pauli corrected) energies from $active.j
    read_eng_blocks(active, "j") do read_block
        [read_block()
         for jj in maximum(j2):-2:minimum(j2)]
    end
end

function atsp_clean()
    # Clean up temporary files that take up a lot of space
    for f in filter(f -> contains(f, ".lst"), readdir())
        rm(f)
    end
end

export csfgenerate_input, csfgenerate, nonh, hf, mchf, atsp_save_run, breit_pauli,
atsp_cp_wfn, read_hf_eng, read_mchf_eng, read_breit_pauli_eng, atsp_clean
