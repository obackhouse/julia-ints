# Compute the 1e overlap integrals

include("gto.jl")
include("expansion.jl")


function _int1e_ovlp_integral(ab::ContractedGaussianPair)
    # Overlap integral between two contracted GTOs
    
    s = 0.0
    la, ma, na = ab.a.lmn
    lb, mb, nb = ab.b.lmn

    for i = 1:ab.a.size
        for j = 1:ab.b.size
            sx = ab.EABx[1,i,j]
            sy = ab.EABy[1,i,j]
            sz = ab.EABz[1,i,j]

            s += ab.D[i,j] * ab.q[i,j]^1.5 * sx * sy * sz
        end
    end

    s * (2.0 * π)^1.5
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
