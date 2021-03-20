# Compute the Hermite polynomial coefficients

using EllipsisNotation

# PC::Vector{Float64}     # P-C coordinates
# FnT::AbstractVector     # Boys function values with (-2p)^n factor


function base(t, u, v, n, PC, FnT)
    #  n             n
    # R      = (-2 p)  F (T)
    #  0,0,0            n

    Rtuv = FnT[.., n+1]

    Rtuv
end


function oob(t, u, v, n, PC, FnT)
    #  n       
    # R      = 0.0  if t<0 or u<0 or v<0
    #  t,u,v   

    Rtuv = zeros(size(PC)[1], size(PC)[2])

    Rtuv
end


function vertical_t(t, u, v, n, PC, FnT)
    #  n              n+1            n+1
    # R      = (t-1) R        + X   R
    #  t,u,v          t-2,u,v    PC  t-1,u,v

    Rtuv = hermite(t-2, u, v, n+1, PC, FnT) * (t-1) .+
           hermite(t-1, u, v, n+1, PC, FnT) .* view(PC, .., 1)

    Rtuv
end


function vertical_u(t, u, v, n, PC, FnT)
    #  n              n+1            n+1
    # R      = (u-1) R        + Y   R
    #  t,u,v          t,u-2,v    PC  t,u-1,v

    Rtuv = hermite(t, u-2, v, n+1, PC, FnT) * (u-1) .+
           hermite(t, u-1, v, n+1, PC, FnT) .* view(PC, .., 2)

    Rtuv
end


function vertical_v(t, u, v, n, PC, FnT)
    #  n              n+1            n+1
    # R      = (v-1) R        + Z   R
    #  t,u,v          t,u,v-2    PC  t,u,v-1
    
    Rtuv = hermite(t, u, v-2, n+1, PC, FnT) * (v-1) .+
           hermite(t, u, v-1, n+1, PC, FnT) .* view(PC, .., 3)

    Rtuv
end


function hermite(t, u, v, n, PC, FnT)
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
