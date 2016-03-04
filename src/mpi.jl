global mpi = false
global mpi_np = 1
global mpi_tmp = ""

function init_mpi(np, tmp)
    global mpi = true
    global mpi_np = np
    global mpi_tmp = tmp
    ENV["MPI_TMP"] = mpi_tmp
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

export init_mpi, clean_mpi_tmp, mpi_cmd
