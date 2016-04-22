using Graphs
using AtomicLevels

type LevelVertex
    config::Config
    term::Term
    calc_props::Dict
    exp_eng
    wfn
end
import Graphs.attributes
function attributes(l::LevelVertex, g::AbstractGraph)
    Dict{UTF8String,Any}("label" => "$(l.config) $(l.term)")
end

LevelVertex(config::Config, term::Term,
      exp_eng = Inf;
      kwargs...) = LevelVertex(config, term, Dict(kwargs), exp_eng, nothing)

level_graph_init = () -> graph(Vector{LevelVertex}(),Vector{Edge{LevelVertex}}())
level_graph = level_graph_init()

import Base.print, Base.show
function print(io::IO, lg::typeof(level_graph_init()))
    map(reverse(topological_sort_by_dfs(level_graph))) do l
        println(l.config, " ", l.term)
    end
end


function add_level!(dep::Union{LevelVertex,Void}, args...; kwargs...)
    s = LevelVertex(args...; kwargs...)
    global level_graph
    add_vertex!(level_graph, s)
    dep != nothing && add_edge!(level_graph, s, dep)
    s
end
add_level!(args...; kwargs...) = add_level!(nothing, args...; kwargs...)

function add_levels!(dep::Union{LevelVertex,Void}, config::Config, args...; kwargs...)
    ts = terms(config)
    for i in eachindex(ts)
        a = copy(args)
        if length(a) > 0 && length(a[1]) > 1
            a = (a[1][i],a[2:end]...)
        end
        dep = add_level!(dep, config, ts[i], a...; kwargs...)
    end
    dep
end
add_levels!(args...; kwargs...) = add_levels!(nothing, args...; kwargs...)

function clear_levels!()
    global level_graph
    level_graph = level_graph_init()
end

# Used to perform calculation using MPI codes if that is
# requested. However, certain level calculations may be explicitly
# marked as unsuitable for MPI codes (e.g. gst of He II), so these are
# masked.
function maybe_mpi_run(g::Function, d::Dict)
    if get(d, :nompi, false)
        delete!(d, :nompi)
        delete!(d, :mpi_np)
        g(d)
    elseif (np = get(d, :mpi_np, 1)) > 1
        delete!(d, :mpi_np)
        init_mpi(np, joinpath(homedir(), "tmp_mpi")) do
            g(d)
        end
    else
        g(d)
    end
end

function get_and_delete!(d::Dict, sym)
    v = d[sym]
    delete!(d, sym)
    v
end

function get_and_delete!(d::Dict, args...)
    map(args) do a
        get_and_delete!(d, a)
    end
end

get_wfn(wfn::Union{AbstractString,Void}) = wfn
get_wfn(wfn::Vector) = first(wfn)

function calc_level(f::Function, level::LevelVertex, wfn; kwargs...)
    wfn = maybe_mpi_run((d) -> f(level, d, get_wfn(wfn)),
                        merge(Dict(kwargs), level.calc_props))
    println(wfn)
    assert(isfile(wfn))
    wfn
end

# Loop over all levels in topological order such that the output of
# one calculation can be used as input for the next, that has
# requested that.
function calc_levels(f::Function; kwargs...)
    global level_graph

    map(reverse(topological_sort_by_dfs(level_graph))) do l
        e = out_edges(l, level_graph)
        if length(e) > 0
            wfn = first(e).target.wfn
        else
            wfn = nothing
        end

        l.wfn = calc_level(f, l, wfn; kwargs...)
    end
end

function calc_levels_atsp(Z; kwargs...)
    calc_levels(;kwargs...) do l,d,wfn0
        ncorr = get_and_delete!(d, :ncorr)
        hf_mchf_bp(l.config, l.term, Z, ncorr, wfn0; overwrite = false, d...)[2]
    end
end

function plot_level_graph(filename)
    stdin, proc = open(`neato -T$(splitext(filename)[2][2:end]) -o "$filename"`, "w")
    to_dot(AST.level_graph, stdin)
    close(stdin)
    filename
end

export add_level!, add_levels!, clear_levels!, calc_levels, calc_levels_atsp, plot_level_graph
