using AtomicLevels
using Lumberjack

function csfgenerate_input(active::Config,
                           term,
                           nexc,
                           configurations...)
    configurations = map(configurations) do c
        orbital_string(c)
    end
    """$(join(configurations, '\n')) ! Configurations

$(orbital_string(active, false, false, ","))
$(string(term,false)) ! Resulting term
$nexc ! Number of excitations"""
end

function load_csfs(filename::AbstractString)
    tp = r"([0-9]+)([A-Z])([0-9]*)"

    blocks = map(strip, split(readall(filename), "*"))

    csf_blocks = map(blocks[1:end-1]) do block
        block = split(block, "\n")
        closed_orbitals = if block[1][end] != ')'
            orbs = join(["$(o)c" for o in split(block[1])], " ")
            block = block[2:end]
            fill(ref_set_list(orbs))
        else
            Config()
        end

        cfgs = map(block[1:2:end]) do cfg
            [closed_orbitals
             ref_set_list(replace(replace(replace(cfg, "( ", "("), "(", ""), ")", ""))]
        end
        couplings = map(block[2:2:end]) do t
            ms = map(tt -> match(tp, tt), split(t))
            map(ms) do m
                S = (parse(Int, m[1]) - 1)//2
                L = searchindex(ells, lowercase(m[2])) - 1
                sen = length(m[3]) > 0 ? parse(Int, m[3]) : "*"
                L,S,sen
            end
        end
        collect(zip(cfgs,couplings))
    end
end

function write_csfs(io::IO, csfs)
    for block in csfs
        length(block) == 0 && continue
        closed_orbitals = closed(block[1][1])
        write(io, "\n")
        write(io, join([@sprintf("% 3d%s", o[1], ells[o[2]+1])
                        for o in closed_orbitals], " "), "\n")
        for cfg in block
            (closed(cfg[1]) == closed_orbitals) || error("$(cfg[1]) does not have the same orbitals closed as the rest of the block; $(closed_orbitals)")
            for o in open(cfg[1])
                write(io, @sprintf("% 3d%s(% 2d)", o[1], ells[o[2]+1], o[3]))
            end
            write(io, "\n")
            for c in cfg[2]
                write(io, @sprintf("% 2d%s%s",
                                   round(Int, 2c[2]+1),
                                   uppercase(ells[c[1]+1]),
                                   c[3] != "*" ? string(c[3]) : " "))
            end
            write(io, "\n")
        end
        write(io, "*\n")
    end
end

function write_csfs(filename::AbstractString, csfs)
    open(filename, "w") do file
        write_csfs(file, csfs)
    end
end

function count_csfs()
    csfs = split(readall("clist.out"), "\n")
    c_pat = r"([0-9]+)([a-z])\([ 0-9]+\)"
    count(l -> ismatch(c_pat, l), csfs)
end

function filter_csfs(f::Function)
    csfs = map(load_csfs("clist.out")) do block
        filter(block) do cfg
            f(cfg...)
        end
    end
    write_csfs("clist.filtered.out", csfs)
    cpf("clist.out", "clist.orig")
    cpf("clist.filtered.out", "clist.out")
end

function lsgen()
    lsgen_inp = split(readall(open("excitationdata.sh")), "\n")
    lsgen_start = "lsgen << EOF"
    lsgen_end = "EOF"
    i = 1
    j = 1
    while true
        i = findnext(lsgen_inp, lsgen_start, j)
        i == 0 && break
        j = findnext(lsgen_inp, lsgen_end, i)
        pipe_file_run("$atsp/lsgen",
                      join(lsgen_inp[i+1:j-1], "\n"))
        cpf("clist.out", "clist.inp")
    end
end

function csfgenerate(active::Config, csf_filter::Function, lists...)
    info("CSF list generation, active set: $(active)")
    # If core not in list, the value is zero (no core)
    lists = join(map(l -> csfgenerate_input(active, l...), lists), "\ny\n")
    pipe_file_run("$atsp/csfexcitation",
                  """$lists
n ! No more lists

""")
    cpf("csfexcitation.log", "$(string(active)).exc")
    lsgen()

    info("Before filtering: $(count_csfs()) CSFs")
    filter_csfs(csf_filter)
    info("After filtering: $(count_csfs()) CSFs")

    cpf("clist.out", "cfg.inp")

    # Filter out active orbitals not present in the CSF list here
    cfgs = readall("clist.out")
    filter(active) do a
        contains(cfgs, "$(principal_qn(a[1]))$(ells[a[2]+1])")
    end
end

export csfgenerate_input, csfgenerate, load_csfs, write_csfs
