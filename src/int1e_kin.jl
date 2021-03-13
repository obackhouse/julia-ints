# Compute the 1e kinetic energy integrals

include("gto.jl")
include("expansion.jl")


function _int1e_kin_integral(ab::ContractedGaussianPair)
    # Kinetic energy integral between two contracted GTOs

    la, ma, na = ab.a.lmn
    lb, mb, nb = ab.b.lmn
    Lb = lb + mb + nb

    @views begin
        sx = expansion(
                la, lb, 0, 
                ab.KAB[1,:,:], 
                ab.PA[1,:,:], ab.PB[1,:,:], 
                ab.p[:,:], ab.q[:,:]
        )
        sy = expansion(
                ma, mb, 0, 
                ab.KAB[2,:,:], 
                ab.PA[2,:,:], ab.PB[2,:,:], 
                ab.p[:,:], ab.q[:,:]
        )
        sz = expansion(
                na, nb, 0, 
                ab.KAB[3,:,:], 
                ab.PA[3,:,:], ab.PB[3,:,:], 
                ab.p[:,:], ab.q[:,:]
        )

        tx = expansion(
                la, lb+2, 0, 
                ab.KAB[1,:,:], 
                ab.PA[1,:,:], ab.PB[1,:,:], 
                ab.p[:,:], ab.q[:,:]
        )
        ty = expansion(
                ma, mb+2, 0, 
                ab.KAB[2,:,:], 
                ab.PA[2,:,:], ab.PB[2,:,:], 
                ab.p[:,:], ab.q[:,:]
        )
        tz = expansion(
                na, nb+2, 0, 
                ab.KAB[3,:,:], 
                ab.PA[3,:,:], ab.PB[3,:,:], 
                ab.p[:,:], ab.q[:,:]
        )

        if lb > 1
            dx = expansion(
                    la, lb-2, 0, 
                    ab.KAB[1,:,:], 
                    ab.PA[1,:,:], ab.PB[1,:,:], 
                    ab.p[:,:], ab.q[:,:]
            )
        end
        if mb > 1
            dy = expansion(
                    ma, mb-2, 0, 
                    ab.KAB[2,:,:], 
                    ab.PA[2,:,:], ab.PB[2,:,:], 
                    ab.p[:,:], ab.q[:,:]
            )
        end
        if nb > 1
            dz = expansion(
                    na, nb-2, 0, 
                    ab.KAB[3,:,:], 
                    ab.PA[3,:,:], ab.PB[3,:,:], 
                    ab.p[:,:], ab.q[:,:]
            )
        end
    end

    txyz   = ab.β' .* (2 * Lb + 3) .* sx .* sy .* sz
    txyz .-= 2.0 * ab.β'.^2 .* (tx.*sy.*sz .+ sx.*ty.*sz .+ sx.*sy.*tz)

    if lb > 1
        txyz .-= 0.5 * lb * (lb+1) * dx .* sy .* sz
    end
    if mb > 1
        txyz .-= 0.5 * mb * (mb+1) * sx .* dy .* sz
    end
    if nb > 1
        txyz .-= 0.5 * nb * (nb+1) * sx .* sy .* dz
    end

    t = ab.D .* ab.q.^1.5 .* txyz
    t = (2.0 * π)^1.5 * sum(t)

    t
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
