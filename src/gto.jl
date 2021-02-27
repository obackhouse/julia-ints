# Basis functions


struct Primitive
    α::Float64          # exponent
    c::Float64          # coefficient
    lmn::Vector{Int32}  # angular momenta
    A::Vector{Float64}  # coords
    N::Float64          # normalisation
end


struct Contracted
    primitives::Vector{Primitive}
end 


function contract(primitives::AbstractVector{Primitive})
    contracted = Contracted(primitives)

    contracted
end



