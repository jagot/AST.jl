function rnucleus(Z, mass_number, M,
                  spin,
                  dip_moment, quad_moment)
    pipe_file_run("$grasp/rnucleus",
                  """$Z ! Nuclear charge
$mass_number ! Nuclear mass number
n ! Don't revise dimensions
$M ! Mass of neutral atom
$spin ! Nuclear spin [ℏ]
$dip_moment ! Nuclear dipole moment [nuclear magnetons]
$quad_moment ! Nuclear quadruple moment [barn]
n ! Don't revise grid
 
""")
end

function rcsfgenerate_input(active,
                            configurations,
                            j::Range,
                            nexc)
    """$(join(configurations, '\n')) ! Configurations

$active
$(j[1]),$(j[end]) ! jmin,jmax
$nexc ! Number of excitations"""
end

rcsfgenerate_input(active,
                   configurations,
                   term::AbstractString,
                   nexc) = rcsfgenerate_input(active,
                                              configurations,
                                              term_to_2j_range(term),
                                              nexc)
function rcsfgenerate(core, active, lists)
    # If core not in list, the value is zero (no core)
    core = findfirst([:He :Ne :Ar :Kr :Xe :Rn], core)
    lists = join(map(l -> rcsfgenerate_input(active, l...), lists), "\ny\n")
    pipe_file_run("$grasp/rcsfgenerate",
                  """* ! Default ordering
$core ! Core
$lists
n ! No more lists

""")
    cpf("rcsfgenerate.log", "$(active_file(active)).exc")
    cpf("rcsf.out", "rcsf.inp")
end

function rangular()
    pipe_file_run("$grasp/rangular",
                  "y ! Default parameters\n")
end

function rwfnestimate()
    # For interactive runs
    run(`$grasp/rwfnestimate`)
end

function rwfnestimate(estimate...)
    pipe_file_run("$grasp/rwfnestimate",
                  """y ! Default parameters
$(join(map(e -> join(e, '\n'), estimate), "\n*\n"))
*
""")
end

function rmcdhf(vary, spectroscopic, cycles=100)
    no_blocks = length(split(readall("rcsf.inp"), "*"))
    pipe_file_run("$grasp/rmcdhf",
                  """y ! Default parameters
$(join(["1" for i = 1:no_blocks], '\n')) $(no_blocks > 1 ? "\n5" : "")
$vary
$spectroscopic
$cycles
""")
end

function rmcdhf()
    # For interactive runs
    run(`$grasp/rmcdhf`)
end

function grasp_save_run(active)
    run(`$grasp/rsave $(active_file(active))`)
end

function grasp_cp_wfn(a,b)
    mkpath("$b")
    cpf("$a/rwfn.out", "$b/rwfn.out")
end

function parse_eng(level, J, parity,
                   Hartrees,
                   Kaysers,
                   eV)
    Dict(:level => parse(Int, level),
         :J => parse(Int, split(J, "/")[1])//2,
         :parity => parse(Int, "$(parity)1"),
         :Hartrees => parse(Float64, replace(Hartrees, "D", "e")),
         :Kayers => parse(Float64, replace(Kaysers, "D", "e")),
         :eV => parse(Float64, replace(eV, "D", "e")))
end

function read_eng(active)
    open("$(active_file(active)).sum") do file
        for line in enumerate(eachline(file))
            strip(line[2]) == "Eigenenergies:" && break
        end
        readline(file)
        readline(file)
        readline(file)
        energies = []
        for line in eachline(file)
            strip(line) == "" && break
            push!(energies, parse_eng(split(strip(line))...))
        end
        energies
    end
end

function grasp_clean()
    # Clean up temporary files that take up a lot of space
    for f in filter(f -> contains(f, "mcp."), readdir())
        rm(f)
    end
end

export rnucleus, rcsfgenerate_input, rcsfgenerate, rangular, rwfnestimate, rmcdhf, grasp_save_run, grasp_cp_wfn, parse_eng, read_eng, grasp_clean
