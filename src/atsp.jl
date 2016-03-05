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

function atsp_save_run(active)
    cpf("wfn.out", "wfn.inp")
    active = active_file(active)
    cpf("wfn.out", "$active.w")
    # cpf("cfg.out", "$active.c")
    cpf("summry", "$active.s")
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

function read_mchf_eng(active)
    # Read the total (non-relativistic) energy from the first .l file
    # in the directory
    open(filter(f -> contains(f, ".l"), readdir())[1]) do file
        for i = 1:5
            readline(file)
        end
        float(split(readline(file))[2])
    end
end

function atsp_clean()
    # Clean up temporary files that take up a lot of space
    for f in filter(f -> contains(f, ".lst"), readdir())
        rm(f)
    end
end

export csfgenerate_input, csfgenerate, nonh, hf, mchf, atsp_save_run, atsp_cp_wfn, read_hf_eng, read_mchf_eng, atsp_clean
