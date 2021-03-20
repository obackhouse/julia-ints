# Compute the 1e kinetic energy integrals

include("gto.jl")
include("expansion.jl")


function _int1e_kin_integral(ab::ContractedGaussianPair)
    # Kinetic energy integral between two contracted GTOs

    t = 0.0
    la, ma, na = ab.a.lmn
    lb, mb, nb = ab.b.lmn
    Lb = lb + mb + nb

    for i = 1:ab.a.size
        for j = 1:ab.b.size
            sx = expansion(
                    la, lb, 0, 
                    ab.KAB[i,j,1], 
                    ab.PA[i,j,1], 
                    ab.PB[i,j,1], 
                    ab.p[i,j], 
                    ab.q[i,j],
            )
            sy = expansion(
                    ma, mb, 0, 
                    ab.KAB[i,j,2], 
                    ab.PA[i,j,2], 
                    ab.PB[i,j,2], 
                    ab.p[i,j], 
                    ab.q[i,j],
            )
            sz = expansion(
                    na, nb, 0, 
                    ab.KAB[i,j,3], 
                    ab.PA[i,j,3], 
                    ab.PB[i,j,3], 
                    ab.p[i,j], 
                    ab.q[i,j],
            )

            tx = expansion(
                    la, lb+2, 0, 
                    ab.KAB[i,j,1], 
                    ab.PA[i,j,1],
                    ab.PB[i,j,1], 
                    ab.p[i,j], 
                    ab.q[i,j],
            )
            ty = expansion(
                    ma, mb+2, 0, 
                    ab.KAB[i,j,2], 
                    ab.PA[i,j,2], 
                    ab.PB[i,j,2], 
                    ab.p[i,j], 
                    ab.q[i,j],
            )
            tz = expansion(
                    na, nb+2, 0, 
                    ab.KAB[i,j,3], 
                    ab.PA[i,j,3], 
                    ab.PB[i,j,3], 
                    ab.p[i,j], 
                    ab.q[i,j],
            )

            if lb > 1
                dx = expansion(
                        la, lb-2, 0, 
                        ab.KAB[i,j,1], 
                        ab.PA[i,j,1], 
                        ab.PB[i,j,1], 
                        ab.p[i,j], 
                        ab.q[i,j],
                )
            end
            if mb > 1
                dy = expansion(
                        ma, mb-2, 0, 
                        ab.KAB[i,j,2], 
                        ab.PA[i,j,2], 
                        ab.PB[i,j,2], 
                        ab.p[i,j], 
                        ab.q[i,j],
                )
            end
            if nb > 1
                dz = expansion(
                        na, nb-2, 0, 
                        ab.KAB[i,j,3], 
                        ab.PA[i,j,3], 
                        ab.PB[i,j,3], 
                        ab.p[i,j], 
                        ab.q[i,j],
                )
            end

            txyz  = ab.β[j] * (2 * Lb + 3) * sx * sy * sz
            txyz -= 2.0 * ab.β[j]^2 * (tx*sy*sz + sx*ty*sz + sx*sy*tz)

            if lb > 1
                txyz -= 0.5 * lb * (lb+1) * dx * sy * sz
            end
            if mb > 1
                txyz -= 0.5 * mb * (mb+1) * sx * dy * sz
            end
            if nb > 1
                txyz -= 0.5 * nb * (nb+1) * sx * sy * dz
            end

            t += ab.D[i,j] * ab.q[i,j]^1.5 * txyz
        end
    end

    t * (2.0 * π)^1.5
end


function _int1e_kin_matrix(basis::Basis)
    # Kinetic energy integral between all contracted GTOs in the basis

    t = Array{Float64, 2}(undef, basis.size, basis.size)

    for pair in basis.pairs
        tij = _int1e_kin_integral(pair)
        t[pair.i, pair.j] = tij
        t[pair.j, pair.i] = conj(tij)
    end

    t
end
