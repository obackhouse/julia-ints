# Compute the two-electron integrals
#TODO: vectorising with arrays over all 4 functions is a
#      terrible idea for highly-contracted basis sets :(

using Tullio

include("gto.jl")
include("expansion.jl")
include("boys.jl")
include("hermite.jl")


function _int2e_integral(ab::ContractedGaussianPair, cd::ContractedGaussianPair, with_cache=false)
    # Two-electron integral between four contracted GTOs

    vee = 0.0
    la, ma, na = ab.a.lmn
    lb, mb, nb = ab.b.lmn
    lc, mc, nc = cd.a.lmn
    ld, md, nd = cd.b.lmn
    Lab = la + ma + na + lb + mb + nb
    Lcd = lc + mc + nc + ld + md + nd
    Labcd = Lab + Lcd

    μ2 = ones(Labcd+1)
    FnT = zeros(Labcd+1)

    for i = 1:ab.a.size
        for j = 1:ab.b.size
            for k = 1:cd.a.size
                for l = 1:cd.b.size
                    vee_ij = 0.0

                    p_times_q = ab.p[i,j] * cd.p[k,l]
                    p_plus_q = ab.p[i,j] + cd.p[k,l]
                    PQ = view(ab.P, i, j, :) .- view(cd.P, k, l, :)
                    μ = p_times_q / p_plus_q
                    T = μ * sum(PQ .* PQ)

                    fill!(μ2, 1.0)
                    fill!(FnT, 0.0)

                    for n in 2:Labcd+1
                        μ2[n] = μ2[n-1] * (-2.0 * μ)
                    end

                    boys_array!(Labcd, T, FnT)
                    FnT .= FnT .* μ2

                    if with_cache
                        cache_t = Dict()
                        cache_u = Dict()
                        cache_v = Dict()
                        cache_w = Dict()
                        cache_x = Dict()
                        cache_y = Dict()
                    else
                        cache_t = cache_u = cache_v = cache_w = cache_x = cache_y = nothing
                    end

                    for w = 0:lc+ld
                        Ew = expansion(
                                lc, ld, w,
                                cd.KAB[k,l,1],
                                cd.PA[k,l,1],
                                cd.PB[k,l,1],
                                cd.p[k,l], 
                                cd.q[k,l],
                                cache_w,
                        )

                        for x = 0:mc+md
                            Ex = expansion(
                                    mc, md, x,
                                    cd.KAB[k,l,2],
                                    cd.PA[k,l,2],
                                    cd.PB[k,l,2],
                                    cd.p[k,l], 
                                    cd.q[k,l],
                                    cache_x,
                            )
                            Ewx = Ew * Ex

                            for y = 0:nc+nd
                                Ey = expansion(
                                        nc, nd, y,
                                        cd.KAB[k,l,3],
                                        cd.PA[k,l,3],
                                        cd.PB[k,l,3],
                                        cd.p[k,l], 
                                        cd.q[k,l],
                                        cache_y,
                                )
                                Ewxy = Ewx * Ey * ((-1)^(w+x+y))

                                for t = 0:la+lb
                                    Et = expansion(
                                            la, lb, t,
                                            ab.KAB[i,j,1],
                                            ab.PA[i,j,1],
                                            ab.PB[i,j,1],
                                            ab.p[i,j],
                                            ab.q[i,j],
                                            cache_t,
                                    )
                                    Etwxy = Et * Ewxy

                                    for u = 0:ma+mb
                                        Eu = expansion(
                                                ma, mb, u,
                                                ab.KAB[i,j,2],
                                                ab.PA[i,j,2],
                                                ab.PB[i,j,2],
                                                ab.p[i,j],
                                                ab.q[i,j],
                                                cache_u,
                                        )
                                        Etuwxy = Etwxy * Eu

                                        for v = 0:na+nb
                                            Ev = expansion(
                                                    na, nb, v,
                                                    ab.KAB[i,j,3],
                                                    ab.PA[i,j,3],
                                                    ab.PB[i,j,3],
                                                    ab.p[i,j],
                                                    ab.q[i,j],
                                                    cache_v,
                                            )
                                            Rtuvwxy = hermite(t+w, u+x, v+y, 0, PQ, FnT)
                                            vee_ij += Etuwxy * Ev * Rtuvwxy
                                        end
                                    end
                                end
                            end
                        end
                    end

                    vee += ab.D[i,j] * cd.D[k,l] * vee_ij / (p_times_q * p_plus_q^0.5)
                end
            end
        end
    end
    
    vee * 2.0 * π^2.5
end


function _int2e_matrix(basis::Basis)
    # Two-electron integral between all contracted GTOs in the basis

    v = Array{Float64, 4}(undef, basis.size, basis.size, basis.size, basis.size)

    npair = length(basis.pairs)

    for ij = 1:npair
        pij = basis.pairs[ij]
        i, j = pij.i, pij.j

        for kl = 1:ij
            pkl = basis.pairs[kl]
            k, l = pkl.i, pkl.j

            vijkl = _int2e_integral(pij, pkl)
            v[i, j, k, l] = vijkl
            v[j, i, k, l] = vijkl
            v[i, j, l, k] = vijkl
            v[j, i, l, k] = vijkl
            v[k, l, i, j] = conj(vijkl)
            v[k, l, j, i] = conj(vijkl)
            v[l, k, i, j] = conj(vijkl)
            v[l, k, j, i] = conj(vijkl)
        end
    end

    v
end


# Legacy: vectorising with arrays over all four contractions is
# super slow, probably too much allocation of small arrays...
#function _int2e_integral(ab::ContractedGaussianPair, cd::ContractedGaussianPair)
#    # Two-electron integral between four contracted GTOs
#
#    vee = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size)
#    la, ma, na = ab.a.lmn
#    lb, mb, nb = ab.b.lmn
#    lc, mc, nc = cd.a.lmn
#    ld, md, nd = cd.b.lmn
#    Lab = la + ma + na + lb + mb + nb
#    Lcd = lc + mc + nc + ld + md + nd
#    Labcd = Lab + Lcd
#
#    p_times_q = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size)
#    p_plus_q = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size)
#    μ = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size)
#    μ2 = ones(ab.a.size, ab.b.size, cd.a.size, cd.b.size, Labcd+1)
#    PQ = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size, 3)
#    T = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size)
#    FnT = zeros(ab.a.size, ab.b.size, cd.a.size, cd.b.size, Labcd+1)
#
#    @tullio p_plus_q[i,j,k,l] = ab.p[i,j] + cd.p[k,l]
#    @tullio p_times_q[i,j,k,l] = ab.p[i,j] * cd.p[k,l]
#    @tullio PQ[i,j,k,l,x] = ab.P[i,j,x] - cd.P[k,l,x]
#    μ .= p_times_q ./ p_plus_q
#    @tullio T[i,j,k,l] = μ[i,j,k,l] * PQ[i,j,k,l,x] * PQ[i,j,k,l,x]
#
#    for n in 2:Labcd+1
#        μ2[:,:,:,:,n] = view(μ2, :, :, :, :, n-1) .* (-2.0 * μ)
#    end
#
#    for i = 1:ab.a.size
#        for j = 1:ab.b.size
#            for k = 1:cd.a.size
#                for l = 1:cd.b.size
#                    boys_array!(Labcd, T[i,j,k,l], view(FnT, i, j, k, l, :))
#                end
#            end
#        end
#    end
#
#    FnT .= FnT .* μ2
#
#    T = nothing
#    μ2 = nothing
#
#    @views begin
#        for w = 0:lc+ld
#            Ew = expansion(
#                    lc, ld, w,
#                    cd.KAB[:,:,1],
#                    cd.PA[:,:,1],
#                    cd.PB[:,:,1],
#                    cd.p, cd.q,
#            )
#
#            for x = 0:mc+md
#                Ex = expansion(
#                        mc, md, x,
#                        cd.KAB[:,:,2],
#                        cd.PA[:,:,2],
#                        cd.PB[:,:,2],
#                        cd.p, cd.q,
#                )
#                Ewx = Ew .* Ex
#
#                for y = 0:nc+nd
#                    Ey = expansion(
#                            nc, nd, y,
#                            cd.KAB[:,:,3],
#                            cd.PA[:,:,3],
#                            cd.PB[:,:,3],
#                            cd.p, cd.q,
#                    )
#                    Ewxy = Ewx .* Ey * ((-1)^(w+x+y))
#
#                    for t = 0:la+lb
#                        Et = expansion(
#                                la, lb, t,
#                                ab.KAB[:,:,1],
#                                ab.PA[:,:,1],
#                                ab.PB[:,:,1],
#                                ab.p, ab.q,
#                        )
#
#                        for u = 0:ma+mb
#                            Eu = expansion(
#                                    ma, mb, u,
#                                    ab.KAB[:,:,2],
#                                    ab.PA[:,:,2],
#                                    ab.PB[:,:,2],
#                                    ab.p, ab.q,
#                            )
#                            Etu = Et .* Eu
#
#                            for v = 0:na+nb
#                                Ev = expansion(
#                                        na, nb, v,
#                                        ab.KAB[:,:,3],
#                                        ab.PA[:,:,3],
#                                        ab.PB[:,:,3],
#                                        ab.p, ab.q,
#                                )
#                                Rtuvwxy = hermite(t+w, u+x, v+y, 0, PQ, FnT)
#                                Etuv = Etu .* Ev
#
#                                @tullio vee[i,j,k,l] = vee[i,j,k,l] + Etuv[i,j] * Ewxy[k,l] * Rtuvwxy[i,j,k,l]
#                            end
#                        end
#                    end
#                end
#            end
#        end
#    end
#
#    @tullio vee[i,j,k,l] = ab.D[i,j] * cd.D[k,l] * vee[i,j,k,l] /
#            (p_times_q[i,j,k,l] * p_plus_q[i,j,k,l]^0.5)
#    vee = 2.0 * π^2.5 * sum(vee)
#
#    vee
#end
