global mpi = false
global mpi_np = 1
global mpi_tmp = ""
global mpi_disks = ""

function init_mpi(f::Function, np::Int, tmp::AbstractString)
    global mpi = true
    global mpi_np = np
    global mpi_tmp = "$tmp/$(DateTime(now()))"
    ENV["MPI_TMP"] = mpi_tmp
    mkdir(mpi_tmp)
    dirs = []
    for i = 1:np
        mkdir(@sprintf("%s/%03d", mpi_tmp, i-1))
        push!(dirs, "'$mpi_tmp'")
    end

    global mpi_disks = join(dirs, '\n')

    res = f()

    deinit_mpi()

    res
end

function deinit_mpi()
    isdir(mpi_tmp) && rm(mpi_tmp, recursive=true)
    mpi = false
    mpi_np = 1
end

function write_mpi_disks()
    open("disks", "w") do file
        write(file, "'$(pwd())'\n")
        write(file, mpi_disks)
    end
end

function clean_mpi_tmp(patterns...)
    mpi || return
    files = readdir(mpi_tmp)
    for p in patterns
        for f in filter(f -> ismatch(p,f), files)
            rm("$mpi_tmp/$f")
        end
    end
end

function mpi_cmd(cmd::AbstractString)
    # If MPI is enabled, return command for running MPI variant of
    # exe, else, return normal exe
    cmd = (ismatch(r"_mpi$", cmd) ?
           (mpi ?
            `mpirun -np $mpi_np $cmd`
            : cmd[1:end-4])
           : cmd)
    `$cmd`
end

function mpi_run(cmd::Union{Cmd,Base.CmdRedirect})
    mpi && write_mpi_disks()
    run(cmd)
end

export init_mpi, write_mpi_disks, clean_mpi_tmp, mpi_cmd, mpi_run
