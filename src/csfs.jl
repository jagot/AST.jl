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

function count_csfs()
    csfs = split(readall("clist.out"), "\n")
    c_pat = r"([0-9]+)([a-z])\([ 0-9]+\)"
    count(l -> ismatch(c_pat, l), csfs)
end

function filter_csfs(f::Function)
    p = r"([0-9]+)([a-z])\(([ 0-9]+?)\)"
    tp = r"([0-9]+)([A-Z])([0-9]*)"
    open("clist.out") do infile
        open("clist.filtered.out", "w") do outfile
            while !eof(infile)
                line = readline(infile)
                if ismatch(p,line)
                    c = map(matchall(p,line)) do m
                        m = match(p,m)
                        "$(m[1])$(m[2])$(strip(m[3]))"
                    end
                    config = ref_set_list(join(c, " "))
                    term_line = readline(infile)
                    m = match(tp, split(term_line)[end])
                    term = Term(searchindex(ells,lowercase(m[2]))-1,
                                (parse(Int, m[1])-1)//2, parity(config))
                    if !f(config, term)
                        continue
                    end
                    write(outfile, line)
                    write(outfile, term_line)
                else
                    write(outfile, line)
                end
            end
        end
    end
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
