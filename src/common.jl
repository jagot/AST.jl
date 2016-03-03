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
    pipe_file = "$(cmd).inp"
    stdout_file = "$(cmd).stdout"
    open(pipe_file, "w") do file
        write(file, contents)
    end
    run(pipeline(`$cmd`,
                 stdin=pipe_file))
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
