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
    lab, mab, nab = la+lb, ma+mb, na+nb
    Lab = la + ma + na + lb + mb + nb
    Rtuv = Array{Float64, 3}(undef, lab+1, mab+1, nab+1)

    for i = 1:ab.a.size
        for j = 1:ab.b.size
            vne_ij = 0.0
            p2 = ones(Lab+1)

            for n in 2:Lab+1
                p2[n] = p2[n-1] * (-2.0 * ab.p[i,j])
            end
            
            for k = 1:length(atoms)
                FnT = zeros(Lab+1)
                PC = view(ab.P, i, j, :) - atoms[k][2]
                T = ab.p[i,j] * sum(PC .* PC)

                boys_array!(Lab, T, FnT)
                FnT .= FnT .* p2

                Rtuv = populate_hermite(
                        lab, mab, nab,
                        PC, FnT, Rtuv,
                )

                for t = 0:lab
                    for u = 0:mab
                        for v = 0:nab
                            vne_ij += charges[k] *
                                      ab.EABx[t+1,i,j] *
                                      ab.EABy[u+1,i,j] *
                                      ab.EABz[v+1,i,j] *
                                      Rtuv[t+1,u+1,v+1]
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
