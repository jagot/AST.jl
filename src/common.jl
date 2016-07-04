using JLD
using AtomicLevels
using Lumberjack

function dir_run(f::Function, wd::AbstractString)
    Lumberjack.info("Entering $(abspath(wd))")
    mkpath(wd)
    here = pwd()
    cd(wd)
    res = 0
    try
        res = f()
    finally
        cd(here)
    end
    Lumberjack.info("Leaving $(abspath(wd))")
    res
end

function pipe_file_run(cmd::AbstractString,
                       contents::AbstractString)
    println(cmd)
    base_cmd = basename(cmd)
    cmd = mpi_cmd(cmd)

    inp_file = basename("$(base_cmd).inp")
    stdout_file = basename("$(base_cmd).stdout")

    open(inp_file, "w") do file
        write(file, contents)
    end

    info("Executing $(cmd) with input file $(abspath(inp_file))")
    debug("Input:\n$(contents)")
    mpi_run(pipeline(inp_file, cmd, stdout_file))
    info("Finished $(cmd)")
end

cpf(src,dst) = cp(src, dst, remove_destination=true)

y_or_n(tf) = tf ? "y" : "n"

function save_eng(energies)
    jldopen("energies.jld", "w") do file
        write(file, "energies", energies)
    end
end

function load_eng(directory)
    jldopen("$directory/energies.jld") do file
        read(file, "energies")
    end
end

function principal_qn(n)
    if n < 10
        n
    else
        ":;<=>?"[n-9]
    end
end

function orbital_string(c::Config,
                        occ = true, status = true,
                        sep = "";
                        escape_principal_qn = false)
    c = map(c) do orb
        s = "$(escape_principal_qn ? principal_qn(orb[1]) : orb[1])$(ells[orb[2]+1])"
        if occ
            if status
                s = string(s, "($(orb[3]),$(orb[4]))")
            else
                s = string(s, "($(orb[3]))")
            end
        else
            s
        end
    end
    join(c, sep)
end

active_file(active::Config, term::Term) = string(join(map(string, active), "_"), "_", string(term))

uncorrelated(filename) = "$(dirname(dirname(filename)))/0/$(basename(filename))"

export dir_run, pipe_file_run, cpf, y_or_n, save_eng, load_eng, uncorrelated
