using AtomicLevels

# Only capable of providing active sets for configurations with no
# holes below in any ell channel.
function active_set(ref_set::Config, i)
    open_set = filter(orb -> orb[4] == "*", ref_set)
    if i == 0
        return open_set
    end

    # Find all occupied n quantum numbers for each ell
    occ = Dict{Int,Vector{Int}}()
    for o in open_set
        ell = o[2]
        occ[ell] = push!(get(occ, ell, Vector{Int}([])), o[1])
    end

    ell_occ = keys(occ)
    ell_max = maximum(ell_occ)

    trunc_neg = v -> v > 0 ? v : 0

    active = Dict{Char,Int}()
    c = map(0:(ell_max+i)) do ell
        max_occ = maximum(get(occ, ell, [0]))
        # Amount of extra to add, if this ell channel is not occupied at all
        extra = trunc_neg(round(Int, ell>ell_max)*(i-(ell-ell_max)))
        n = max(max_occ + i, ell + extra + 1)
        Orbital((n,ell,1,"*"))
    end
end
