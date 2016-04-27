using Graphs
using Jagot.plotting

function result_table()
    s = map(reverse(topological_sort_by_dfs(level_graph))) do l
        Dict(:heading => "$(latex(l.config)) $(latex(l.term))",
             :config => l.config,
             :term => l.term,
             :energies => load_eng(dirname(dirname(l.wfn))),
             :exp_eng => l.exp_eng)
    end

    sort!(s, by=v -> minimum(v[:energies]))

    ml = maximum([size(v[:energies],2) for v in s])

    function extend(v,ml)
        m = size(v,1)
        if m == ml
            v
        else
            [v[1:end-1,:]
             repmat(v[end,:],ml-m+1,1)]
        end
    end

    latex_number(::Void) = ""
    latex_number(n::Number) = "\\($(latex_base10(n))\\)"

    function latex_numbers(v::Vector)
        map(v) do e
            latex_number(e)
        end
    end

    function latex_numbers(m::Matrix)
        map(m) do v
            latex_number(v)
        end
    end

    labels = hcat([[v[:heading] repmat([""], 1, size(v[:energies],1)-1)] for v in s]...)
    labels2 = map(s) do v
        ["(MC)HF" map(latex, levels(v[:config], v[:term]))...]
    end
    labels2 = hcat(labels2...)

    # Breit--Pauli yields the energies in order of decreasing J.
    eng = map(s) do v
        ee = v[:energies]'
        extend([ee[:,1] ee[:,end:-1:2]], ml)
    end
    eng = hcat(eng...)
    deng = copy(eng)
    deng[2:end,:] -= deng[1:end-1,:]
    deng2 = broadcast(-, eng, eng[:,indmin(eng[end,:])])

    exp_eng = map(s) do v
        ee = v[:exp_eng][1]
        if length(ee) == 1
            [ee repmat([nothing], 1, size(v[:energies],1)-1)]
        else
            [nothing ee...]
        end
    end
    exp_eng = hcat(exp_eng...)

    layers = 0:(size(eng,1)-1)

    spacer = repmat([""],1,size(labels,2)+1)

    ["Layer" labels
     "" labels2
     spacer
     layers latex_numbers(eng)
     spacer
     layers latex_numbers(deng)
     spacer
     layers latex_numbers(deng2)
     spacer
     "Exp" latex_numbers(exp_eng)]
end

export result_table
