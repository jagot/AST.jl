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
    abs(2L-S2):(2L+S2)
end

export ells, degeneracy, split_ref_set, term_to_2j_range
