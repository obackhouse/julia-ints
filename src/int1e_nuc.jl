# Compute the nuclear repulsion integrals

using LinearAlgebra
using Tullio

include("gto.jl")
include("expansion.jl")
include("boys.jl")
include("hermite.jl")


function _int1e_nuc_integral(ab::ContractedGaussianPair, atoms::Vector{Tuple{String, Vector{Float64}}}, charges::Vector{Float64})
    # Nuclear repulsion integral between two contracted GTOs

    vne = 0.0
    la, ma, na = ab.a.lmn
    lb, mb, nb = ab.b.lmn
    Lab = la + ma + na + lb + mb + nb

    mtp = ones(Lab+1, ab.a.size, ab.b.size)
    vne = zeros(ab.a.size, ab.b.size)

    for n in 2:Lab+1
        mtp[n,:,:] = view(mtp, n-1, :, :) .* (-2.0 * ab.p)
    end
    
    for k = 1:length(atoms)
        coords = atoms[k][2]
        FnT = zeros(Lab+1, ab.a.size, ab.b.size)
        PC = zeros(3, ab.a.size, ab.b.size)
        T = zeros(ab.a.size, ab.b.size)

        @tullio PC[x,i,j] = ab.P[x,i,j] - coords[x]
        @tullio T[i,j] = ab.p[i,j] * PC[x,i,j] * PC[x,i,j]

        for i = 1:ab.a.size
            for j = 1:ab.b.size
                boys_array!(Lab, T[i,j], view(FnT, :, i, j))
            end
        end

        FnT .= FnT .* mtp

        @views begin
            for t = 0:la+lb
                Et = expansion(
                        la, lb, t, 
                        ab.KAB[1,:,:],
                        ab.PA[1,:,:],
                        ab.PB[1,:,:],
                        ab.p, ab.q,
                )

                for u = 0:ma+mb
                    Eu = expansion(
                            ma, mb, u, 
                            ab.KAB[2,:,:],
                            ab.PA[2,:,:],
                            ab.PB[2,:,:],
                            ab.p, ab.q,
                    )

                    for v = 0:na+nb
                        Ev = expansion(
                                na, nb, v, 
                                ab.KAB[3,:,:],
                                ab.PA[3,:,:],
                                ab.PB[3,:,:],
                                ab.p, ab.q,
                        )
                        Rtuv = hermite(t, u, v, 0, PC, FnT)

                        vne .+= charges[k] * Et .* Eu .* Ev .* Rtuv
                    end
                end
            end
        end
    end

    vne = ab.D .* ab.q .* vne
    vne = -4.0 * π * sum(vne)

    vne
end


function _int1e_nuc_matrix(basis::Basis)
    # Nuclear repulsion integral between all contracted GTOs in the basis

    v = Array{Float64, 2}(undef, basis.size, basis.size)

    for pair in basis.pairs
        vij = _int1e_nuc_integral(pair, atoms, basis.charges)
        v[pair.i, pair.j] = vij
        v[pair.j, pair.i] = conj(vij)
    end

    v
end
