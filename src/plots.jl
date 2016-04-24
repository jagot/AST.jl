using InterpNd
using PyPlot

function add_partial_waves(data, orbitals)
    ml = maximum(map(d -> size(d,2), data))
    mr = maximum(map(d -> maximum(d[1,:]), data))
    rsq = linspace(0,mr,ml)
    core = zeros(ml)
    for i in eachindex(orbitals)
        w = degeneracy(orbitals[i])
        ii = interp1d(vec(data[i][2,:]))
        core += w*abs(ii(vec(data[i][1,:]),rsq)).^2
    end

    rsq,core
end

function plot_wfn(config::Config, term::Term, ncorr)
    filename = active_file(config, term)
    conf_orbs = map(config) do o
        Orbital((o[1],o[2],1,"*"))
    end
    dir = "$filename/$ncorr"
    sep = "    0.0000       0.000"
    cp = r"([0-9]+)([a-z])\(([ 0-9]+?)\)"
    orbitals,data = dir_run(dir) do
        clist = map(filter(l -> ismatch(cp, l), eachline("$(filename).c"))) do l
            map(matchall(cp,l)) do m
                ref_set_list("$(m[1])$(m[2])")[1]
            end
        end
        orbitals = sort(unique(vcat(conf_orbs...,clist...)))
        rd,wr = redirect_stdin()
        for o in orbitals
            write(wr, " \n")
        end
        run(pipeline(rd, `$(atsp)/plotw $(filename).w`))
        orbitals,map(split(readall("plot.dat"), sep)[2:end]) do pw
            v = map([sep; split(pw, "\n")[2:end-4]]) do l
                map(f -> parse(Float64, f), split(l))
            end
            hcat(v...)
        end
    end
    subplot(211)
    for i in eachindex(orbitals)
        PyPlot.plot(data[i][1,:]',data[i][2,:]'.^2,label=string(orbitals[i]))
    end
    legend(framealpha=0.75)
    margins(0,0.1)
    subplot(212)
    PyPlot.plot(add_partial_waves(data,orbitals)...)
    margins(0,0.1)
end

export plot_wfn
