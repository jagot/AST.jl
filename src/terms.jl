# http://download.springer.com/static/pdf/543/bbm%253A978-1-4419-8182-0%252F1.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Fbook%2Fbbm%3A978-1-4419-8182-0%2F1&token2=exp=1457644468~acl=%2Fstatic%2Fpdf%2F543%2Fbbm%25253A978-1-4419-8182-0%25252F1.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Fbook%252Fbbm%253A978-1-4419-8182-0%252F1*~hmac=4c1a0196497a49ccd369ab63df10c876e3c8bd9ba48e2654fadc300f9531186b

add_electron(J::Void, ell) = ell:ell
add_electron(J::Range, ell) = abs(minimum(J)-ell):(maximum(J)+ell)

function terms(ref_set)

    L = nothing
    S = nothing
    for orb in ref_set
        ell = searchindex(ells,orb[1][end])-1
        g_ell = degeneracy(ell)
        g_ell == orb[2] && continue
        for i = (orb[2] <= g_ell/2 ? (1:orb[2]) : (1:(g_ell-orb[2])))
            L = add_electron(L, ell)
            S = add_electron(S, 1//2)
        end
    end
    if L === nothing
        L = 0:0
        S = 0:0
    end
    vcat([["$(round(Int, ESS2P1))$(uppercase(ells[ELL+1]))" for ESS2P1 in 2S + 1] for ELL in L]...)
end

export terms
