# Compute the two-electron integrals

using Tullio

include("gto.jl")
include("expansion.jl")
include("boys.jl")
include("hermite.jl")


function _int2e_integral_ssss(ab::ContractedGaussianPair, cd::ContractedGaussianPair)
    # Two-electron integral between four contracted s-type GTOs
    
    vee = 0.0

    for i = 1:ab.a.size
        for j = 1:ab.b.size
            for k = 1:cd.a.size
                for l = 1:cd.b.size
                    p_times_q = ab.p[i,j] * cd.p[k,l]
                    p_plus_q = ab.p[i,j] + cd.p[k,l]
                    PQ = view(ab.P, i, j, :) .- view(cd.P, k, l, :)
                    T = (p_times_q / p_plus_q) * sum(PQ .* PQ)
                    FnT = boys(0, T)

                    vee += ab.D[i,j] * cd.D[k,l] * FnT / (p_times_q * p_plus_q^0.5) *
                           ab.KAB[i,j,1] * ab.KAB[i,j,2] * ab.KAB[i,j,3] *
                           cd.KAB[k,l,1] * cd.KAB[k,l,2] * cd.KAB[k,l,3]
                end
            end
        end
    end

    vee * 2.0 * π^2.5
end


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

    if Labcd == 0
        return _int2e_integral_ssss(ab, cd)
    end

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


function _get_metric(basis::Basis, v::AbstractArray{T, N}) where {T, N}
    # Get the Q_ij metric, i.e. Schwarz inequality = √(ij|ij)

    q = Array{Float64, 2}(undef, basis.size, basis.size)

    npair = length(basis.pairs)

    for ij = 1:npair
        pij = basis.pairs[ij]
        i, j = pij.i, pij.j

        vijij = @inbounds _int2e_integral(pij, pij)
        v[i,j,i,j] = vijij
        v[j,i,i,j] = vijij
        v[i,j,j,i] = vijij
        v[j,i,j,i] = vijij

        qij = vijij^0.5
        q[i,j] = qij
        q[j,i] = conj(qij)
    end

    q
end


function _int2e_matrix(basis::Basis, tol=1e-16)
    # Two-electron integral between all contracted GTOs in the basis

    v = Array{Float64, 4}(undef, basis.size, basis.size, basis.size, basis.size)
    q = _get_metric(basis, v)

    npair = length(basis.pairs)

    for ij = 1:npair
        pij = basis.pairs[ij]
        i, j = pij.i, pij.j
        qij = q[i,j]

        for kl = 1:(ij-1)  # ij==kl handled in metric
            pkl = basis.pairs[kl]
            k, l = pkl.i, pkl.j
            qkl = q[k,l]

            if (q[i,j] * q[k,l]) < tol
                continue
            end

            vijkl = @inbounds _int2e_integral(pij, pkl)
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
