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

function j2_to_jstr(j2)
    jstr = i -> round(Int,i) == i ? round(Int, i) : "$(num(i))/$(den(i))"
    s = join(map(jstr, (minimum(j2):2:maximum(j2))//2), ", ")
    "[$s]"
end

export ells, degeneracy, split_ref_set, term_to_2j_range, j2_to_jstr
