# Compute the 1e overlap integrals

include("gto.jl")
include("expansion.jl")


function _int1e_ovlp_integral(ab::ContractedGaussianPair)
    # Overlap integral between two contracted GTOs
    
    la, ma, na = ab.a.lmn
    lb, mb, nb = ab.b.lmn

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
    end

    s = ab.D .* ab.q.^1.5 .* sx .* sy .* sz
    s = (2.0 * π)^1.5 * sum(s)

    s
end


function _int1e_ovlp_matrix(basis::Basis)
    # Overlap integral between all contracted GTOs in the basis

    s = Array{Float64, 2}(undef, basis.size, basis.size)

    for pair in basis.pairs
        sij = _int1e_ovlp_integral(pair)
        s[pair.i, pair.j] = sij
        s[pair.j, pair.i] = conj(sij)
    end

    s
end
