using AtomicLevels

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

function csfgenerate(active::Config, csf_filter::Function, lists...)
    # If core not in list, the value is zero (no core)
    lists = join(map(l -> csfgenerate_input(active, l...), lists), "\ny\n")
    pipe_file_run("$atsp/csfexcitation",
                  """$lists
n ! No more lists

""")
    cpf("csfexcitation.log", "$(string(active)).exc")
    lsgen_inp = split(readall(open("excitationdata.sh")), "\n")
    pipe_file_run("$atsp/lsgen",
                  join(lsgen_inp[3:end-2], "\n"))

    filter_csfs(csf_filter)

    cpf("clist.out", "cfg.inp")

    # Filter out active orbitals not present in the CSF list here
    cfgs = readall("clist.out")
    filter(active) do a
        contains(cfgs, "$(a[1])$(ells[a[2]+1])")
    end
end
