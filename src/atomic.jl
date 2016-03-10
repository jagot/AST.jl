ells = "spdfghiklmnoqrtuvwxyz"
degeneracy(ell) = 2(2ell + 1)

function split_ref_set(ref_set)
    pattern = r"([0-9]+)([a-z])\(([0-9]+),(.)\)"
    matches = matchall(pattern, ref_set)
    orb = m -> (parse(Int, m[1]),
               searchindex(ells, m[2]) - 1,
               parse(Int, m[3]),
               m[4])
    [orb(match(pattern,m)) for m in matches]
end

function term_to_2j_range(term)
    m = match(r"([0-9]+)([A-Z])", term)
    S2 = parse(Int, m[1]) - 1
    L = searchindex(ells, lowercase(m[2])) - 1
    abs(2L-S2):2:(2L+S2)
end

jstr(i) = round(Int,i) == i ? round(Int, i) : "$(num(i))/$(den(i))"

function j2_to_jstr(j2)
    s = join(map(jstr, (minimum(j2):2:maximum(j2))//2), ", ")
    "[$s]"
end

# Only capable of providing active sets for configurations with no
# holes below in any ell channel.
function active_set(ref_set, i)
    open_set = filter(orb -> orb[3] != "c", ref_set)

    occ = Dict{Char,Vector{Int}}()
    for o in open_set
        ell = o[1][end]
        occ[ell] = push!(get(occ, ell, Vector{Int}([])), parse(Int, o[1][1:end-1]))
    end

    ell_ind = Vector{Int}([searchindex(ells,ell) for ell in keys(occ)])
    ell_min,ell_max = minimum(ell_ind),maximum(ell_ind)

    trunc_neg = v -> v > 0 ? v : 0

    active = Dict{Char,Int}()
    c = map(1:(ell_max+i)) do ellp1
        ell = ells[ellp1]
        max_occ = maximum(get(occ, ell, [0]))
        # Amount of extra to add, if this ell channel is not occupied at all
        extra = trunc_neg(round(Int, ellp1>ell_max)*(i-(ellp1-ell_max)))
        n = max(max_occ + i, ellp1 + extra)
        "$n$ell"
    end
    join(c, ",")
end

export ells, degeneracy, split_ref_set, term_to_2j_range, jstr, j2_to_jstr, active_set
