# Basis functions
#FIXME: last element in basis set file is not being parsed

include("utils.jl")
include("ptable.jl")
include("expansion.jl")


ang_keys = ["S", "P", "D", "F", "G", "H", "I", "J", "K"]


abstract type Gaussian end


struct PrimitiveGaussian <: Gaussian
    α::Float64          # exponent
    c::Float64          # coefficient
    lmn::Vector{Int64}  # angular momenta
    A::Vector{Float64}  # coords
    N::Float64          # normalisation
end


struct ContractedGaussian <: Gaussian
    α::Vector{Float64}  # exponent
    c::Vector{Float64}  # coefficient
    lmn::Vector{Int64}  # angular momenta
    A::Vector{Float64}  # coords
    N::Vector{Float64}  # normalisation
    size::Int64         # number of primitives
end


struct ContractedGaussianPair <: Gaussian
    i::Int64                # index of function a in basis
    j::Int64                # index of function b in basis
    a::ContractedGaussian   # contracted gaussian a
    b::ContractedGaussian   # contracted gaussian b
    α::Array{Float64}       # exponent of a
    β::Array{Float64}       # exponent of b
    p::Array{Float64}       # α+β
    q::Array{Float64}       # 1/(2(α+β))
    D::Array{Float64}       # product c and N of a and b
    P::Array{Float64}       # gaussian product centre
    PA::Array{Float64}      # P-A coordinates
    PB::Array{Float64}      # P-B coordinates
    AB::Array{Float64}      # A-B coordinates
    KAB::Array{Float64}     # exponential term
    EABx::Array{Float64}    # expansion coefficients in x direction
    EABy::Array{Float64}    # expansion coefficients in y direction
    EABz::Array{Float64}    # expansion coefficients in z direction
end


struct Basis
    gaussians::Vector{ContractedGaussian}
    pairs::Vector{ContractedGaussianPair}
    atoms::Vector{Tuple{String, Vector}}
    charges::Vector{Float64}
    size::Int64
    name::String
    dict::AbstractDict
end


function contract(primitives::AbstractVector{PrimitiveGaussian})
    # Convert a list of primitive gaussian to a contract gaussian
    
    α = Vector{Float64}() 
    c = Vector{Float64}()
    N = Vector{Float64}()
    lmn = primitives[1].lmn
    A = primitives[1].A
    
    for i = 1:length(primitives)
        push!(α, primitives[i].α)
        push!(c, primitives[i].c)
        push!(N, primitives[i].N)
        
        @assert isapprox(A, primitives[i].A)
        @assert isequal(lmn, primitives[i].lmn)
    end

    contracteds = ContractedGaussian(α, c, lmn, A, N, length(α))

    contracteds
end


function normalise(α::Vector{Float64}, lmn::Vector{Int64})
    # Normalise a primitive gaussian and set the respective const 
    
    l, m, n = lmn
    L = l + m + n
    factor = (2.0 / π)^0.75

    facl = factorial_2(2*l-1)
    facm = factorial_2(2*m-1)
    facn = factorial_2(2*n-1)

    N = factor * 2^L * α.^(0.25 * (2 * L + 3)) * (facl * facm * facn)^-0.5

    N
end


function build_pairs(gaussians::AbstractVector{ContractedGaussian})
    # Build the contracted gaussian pairs

    pairs = Vector{ContractedGaussianPair}()
    ngto = length(gaussians)

    for i = 1:ngto
        a = gaussians[i]
        na = length(a.α)
        for j = 1:i
            b = gaussians[j]
            nb = length(b.α)   

            α = a.α
            β = b.α
            A = a.A
            B = b.A
            p = α .+ β'
            q = 1.0 ./ (2.0 * p)
            D = (a.c .* a.N) .* (b.c .* b.N)'

            P = Array{Float64}(undef, na, nb, 3)
            PA = Array{Float64}(undef, na, nb, 3)
            PB = Array{Float64}(undef, na, nb, 3)
            AB = Array{Float64}(undef, na, nb, 3)
            KAB = Array{Float64}(undef, na, nb, 3)
            EABx = Array{Float64}(undef, a.lmn[1]+b.lmn[1]+1, na, nb)
            EABy = Array{Float64}(undef, a.lmn[2]+b.lmn[2]+1, na, nb)
            EABz = Array{Float64}(undef, a.lmn[3]+b.lmn[3]+1, na, nb)

            @views begin
                for k = 1:na
                    for l = 1:nb
                        P[k,l,:] .= (α[k] * A .+ β[l] * B) ./ p[k,l]
                        PA[k,l,:] .= P[k,l,:] .- A
                        PB[k,l,:] .= P[k,l,:] .- B
                        AB[k,l,:] .= A .- B
                        ζ = α[k] * β[l] / p[k,l]
                        KAB[k,l,:] .= exp.(-ζ .* AB[k,l,:].^2)

                        populate_expansion(
                               a.lmn[1],
                               b.lmn[1],
                               a.lmn[1]+b.lmn[1],
                               KAB[k,l,1],
                               PA[k,l,1],
                               PB[k,l,1],
                               p[k,l],
                               q[k,l],
                               view(EABx, :, k, l),
                        )

                        populate_expansion(
                               a.lmn[2],
                               b.lmn[2],
                               a.lmn[2]+b.lmn[2],
                               KAB[k,l,2],
                               PA[k,l,2],
                               PB[k,l,2],
                               p[k,l],
                               q[k,l],
                               view(EABy, :, k, l),
                        )

                        populate_expansion(
                               a.lmn[3],
                               b.lmn[3],
                               a.lmn[3]+b.lmn[3],
                               KAB[k,l,3],
                               PA[k,l,3],
                               PB[k,l,3],
                               p[k,l],
                               q[k,l],
                               view(EABz, :, k, l),
                        )
                    end
                end
            end

            pair = ContractedGaussianPair(i, j, a, b, α, β, p, q,
                                          D, P, PA, PB, AB, KAB,
                                          EABx, EABy, EABz)
            push!(pairs, pair)
        end
    end

    pairs
end


function get_lmn(ang)
    # Take an angular momentum char i.e. S, P and return a vector
    # of vectors with each possible l,m,n angular momenta

    lmn = Vector{Vector{Int64}}()
    L = findall(x -> x == ang, ang_keys)[1] - 1

    for l = 0:L
        for m = 0:(L-l)
            n = L - l - m
            push!(lmn, [l,m,n])
        end
    end

    lmn
end


function parse_basis(atoms::Vector{Tuple{String, Vector{Float64}}}, basis::String)
    # Parse the basis set file and return a Basis structs, complete
    # with functions and pairs.

    gtos = Vector{ContractedGaussian}()
    pairs = Vector{ContractedGaussianPair}()

    basis_file_name = replace(lowercase(basis), "-" => "")
    basis_file_name = replace(basis_file_name, "_" => "")
    path = String(@__DIR__) * "/basis/" * basis_file_name * ".dat"
    lines = readlines(path)

    deleteat!(lines, findall(x -> lstrip(rstrip(x)) == "", lines))
    deleteat!(lines, findall(x -> startswith(x, "#"), lines))

    #           atom    prims        ang     coeffs           exps
    dict = Dict{String, Vector{Tuple{String, Vector{Float64}, Vector{Float64}}}}()
    prims = Vector{Tuple{String, Vector{Float64}, Vector{Float64}}}()
    c = Vector{Float64}()
    c2 = Vector{Float64}()
    α = Vector{Float64}()
    ang = ""
    atom = ""

    i = 1
    while i < (length(lines)+1)
        line = lines[i]
        if split(line)[1] in ang_keys || split(line)[1] == "L"
            if i != 1
                if ang != "L"
                    push!(prims, (ang, c, α))
                    c = Vector{Float64}()
                    α = Vector{Float64}()
                else
                    push!(prims, ("S", c, α))
                    push!(prims, ("P", c2, α))
                    c = Vector{Float64}()
                    c2 = Vector{Float64}()
                    α = Vector{Float64}()
                end
            end
            ang, atom_next = split(line)

            if i != 1 && atom != atom_next
                push!(dict, atom => prims)
                prims = Vector{Tuple{String, Vector{Float64}, Vector{Float64}}}()
            end
            atom = atom_next
        else
            line = replace(line, "D" => "E")

            if length(split(line)) == 2
                ci, αi = parse.(Float64, split(line))
                push!(c, ci)
                push!(α, αi)
            else
                cs, αi, cp = parse.(Float64, split(line))
                push!(c, cs)
                push!(α, αi)
                push!(c2, cp)
            end
        end
        i += 1
    end

    for (i, (atom, coord)) in enumerate(atoms)
        A = coord
        for (ang, coeffs, exps) in dict[atom]
            for lmn in get_lmn(ang)
                pgtos = Vector{PrimitiveGaussian}()
                norms = normalise(exps, lmn)

                for (c, α, N) in zip(coeffs, exps, norms)
                    pgto = PrimitiveGaussian(α, c, lmn, A, N)
                    push!(pgtos, pgto)
                end

                gto = contract(pgtos)
                push!(gtos, gto)
            end
        end
    end

    charges = [get_proton_number(atom[1]) for atom in atoms]
    pairs = build_pairs(gtos)
    basis = Basis(gtos, pairs, atoms, charges, length(gtos), basis, dict)

    basis
end


function summarise(basis::Basis)
    # Print a summary of the basis

    for i = 1:basis.size
        gto = basis.gaussians[i]
        @show i, gto.c, gto.α, gto.N, gto.A
    end
end
