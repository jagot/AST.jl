using Graphs

type State
    ref_set::Vector
    calc_props::Dict
    exp_eng
    wfn
end

State(ref_set::Vector,
      exp_eng = Inf;
      kwargs...) = State(ref_set, Dict(kwargs), exp_eng, nothing)

state_graph = graph(Vector{State}(),Vector{Edge{State}}())

function add_state!(dep::Union{State,Void}, args...; kwargs...)
    s = State(args...; kwargs...)
    global state_graph
    add_vertex!(state_graph, s)
    dep != nothing && add_edge!(state_graph, s, dep)
    s
end
add_state!(args...; kwargs...) = add_state!(nothing, args...; kwargs...)

# Used to perform calculation using MPI codes if that is
# requested. However, certain state calculations may be explicitly
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

function calc_state(f::Function, state::State, wfn; kwargs...)
    wfn = maybe_mpi_run((d) -> f(state, d, get_wfn(wfn)),
                        merge(Dict(kwargs), state.calc_props))
    assert(all(map(isfile, wfn)))
    wfn
end

# Loop over all states in topological order such that the output of
# one calculation can be used as input for the next, that has
# requested that.
function calc_states(f::Function; kwargs...)    
    global state_graph
    
    map(reverse(topological_sort_by_dfs(state_graph))) do s
        e = out_edges(s, state_graph)
        if length(e) > 0
            wfn = first(e).target.wfn
        else
            wfn = nothing
        end

        s.wfn = calc_state(f, s, wfn; kwargs...)
    end
end

function calc_states_atsp(Z; kwargs...)
    calc_states(;kwargs...) do state,d,wfn0
        ncorr = get_and_delete!(d, :ncorr)
        map(terms(state.ref_set)) do term
            # Assigning to wfn0 ensures that the result of one term is
            # used as initial guess for the next term            
            wfn0 = hf_mchf_bp(state.ref_set, term, Z, ncorr, wfn0; overwrite = false, d...)[2]
        end
    end
end

export add_state!, calc_states, calc_states_atsp
