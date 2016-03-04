using JLD

function dir_run(f::Function, wd::AbstractString)
    mkpath(wd)
    here = pwd()
    cd(wd)
    res = 0
    try
        res = f()
    finally
        cd(here)
    end
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

    run(pipeline(inp_file, cmd, stdout_file))
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

export dir_run, pipe_file_run, cpf, y_or_n, save_eng, load_eng
