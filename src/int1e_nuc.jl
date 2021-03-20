# Compute the nuclear repulsion integrals

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

    for i = 1:ab.a.size
        for j = 1:ab.b.size
            vne_ij = 0.0
            mtp = ones(Lab+1)

            for n in 2:Lab+1
                mtp[n] = mtp[n-1] * (-2.0 * ab.p[i,j])
            end
            
            for k = 1:length(atoms)
                FnT = zeros(Lab+1)
                PC = view(ab.P, i, j, :) - atoms[k][2]
                T = ab.p[i,j] * sum(PC .* PC)

                boys_array!(Lab, T, FnT)
                FnT .= FnT .* mtp

                @views begin
                    for t = 0:la+lb
                        Et = expansion(
                                la, lb, t, 
                                ab.KAB[i,j,1],
                                ab.PA[i,j,1],
                                ab.PB[i,j,1],
                                ab.p[i,j], 
                                ab.q[i,j],
                        )

                        for u = 0:ma+mb
                            Eu = expansion(
                                    ma, mb, u, 
                                    ab.KAB[i,j,2],
                                    ab.PA[i,j,2],
                                    ab.PB[i,j,2],
                                    ab.p[i,j], 
                                    ab.q[i,j],
                            )

                            for v = 0:na+nb
                                Ev = expansion(
                                        na, nb, v, 
                                        ab.KAB[i,j,3],
                                        ab.PA[i,j,3],
                                        ab.PB[i,j,3],
                                        ab.p[i,j], 
                                        ab.q[i,j],
                                )
                                Rtuv = hermite(t, u, v, 0, PC, FnT)

                                vne_ij += charges[k] * Et * Eu * Ev * Rtuv
                            end
                        end
                    end
                end
            end

            vne += ab.D[i,j] * ab.q[i,j] * vne_ij
        end
    end

    vne * -4.0 * π
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
