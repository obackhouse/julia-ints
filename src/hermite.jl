# Compute the Hermite polynomial coefficients

# PC::Vector{Float64}     # P-C coordinates
# FnT::AbstractVector     # Boys function values with (-2p)^n factor


function base(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
)
    #  n             n
    # R      = (-2 p)  F (T)
    #  0,0,0            n

    FnT[n+1]
end


function oob(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
)
    #  n       
    # R      = 0.0  if t<0 or u<0 or v<0
    #  t,u,v   

    zero(FnT[1])
end


function vertical_t(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
)
    #  n              n+1            n+1
    # R      = (t-1) R        + X   R
    #  t,u,v          t-2,u,v    PC  t-1,u,v

    Rtuv = hermite(t-2, u, v, n+1, PC, FnT) * (t-1) +
           hermite(t-1, u, v, n+1, PC, FnT) * PC[1]

    Rtuv
end


function vertical_u(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
)
    #  n              n+1            n+1
    # R      = (u-1) R        + Y   R
    #  t,u,v          t,u-2,v    PC  t,u-1,v

    Rtuv = hermite(t, u-2, v, n+1, PC, FnT) * (u-1) +
           hermite(t, u-1, v, n+1, PC, FnT) * PC[2]

    Rtuv
end


function vertical_v(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
)
    #  n              n+1            n+1
    # R      = (v-1) R        + Z   R
    #  t,u,v          t,u,v-2    PC  t,u,v-1
    
    Rtuv = hermite(t, u, v-2, n+1, PC, FnT) * (v-1) +
           hermite(t, u, v-1, n+1, PC, FnT) * PC[3]

    Rtuv
end


function hermite_mmd(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
)
    # Default path for computing terms.

    if (t < 0) || (u < 0) || (v < 0)
        Rtuv = oob(t, u, v, n, PC, FnT)
    elseif t == u == v == 0
        Rtuv = base(t, u, v, n, PC, FnT)
    elseif t == u == 0
        Rtuv = vertical_v(t, u, v, n, PC, FnT)
    elseif t == 0
        Rtuv = vertical_u(t, u, v, n, PC, FnT)
    else
        Rtuv = vertical_t(t, u, v, n, PC, FnT)
    end

    Rtuv
end


function hermite(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
        cache::Nothing=nothing,
)
    # Driver function for hermite terms

    hermite_mmd(t, u, v, n, PC, FnT)
end


function hermite(
        t::Int64, 
        u::Int64, 
        v::Int64, 
        n::Int64, 
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
        cache::Dict{NTuple{4, Int64}, Float64},
)
    # Driver function for hermite terms with caching

    key = (t, u, v, n)
    if haskey(cache, key)
        return cache[key]
    end

    Rtuv = hermite_mmd(t, u, v, n, PC, FnT)

    push!(cache, key => Rtuv)

    Rtuv
end


function populate_hermite(
        tmax::Int64,
        umax::Int64,
        vmax::Int64,
        PC::Vector{Float64}, 
        FnT::Vector{Float64},
        out::Array{Float64, 3},
)
    # Populate an array with hermite coefficients for n = 0 and
    # t = 0 → tmax, u = 0 → umax and v = 0 → vmax

    if (tmax+umax+vmax) > 3  #TODO: I assume this is good for high ang mom?
        cache = Dict{NTuple{4, Int64}, Float64}()
    else
        cache = nothing
    end

    #TODO more efficient path
    for t = 0:tmax
        for u = 0:umax
            for v = 0:vmax
                out[t+1,u+1,v+1] = hermite(t, u, v, 0, PC, FnT, cache)
            end
        end
    end

    out
end
