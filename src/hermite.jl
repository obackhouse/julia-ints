# Compute the Hermite polynomial coefficients


struct HermitePair
    PC::Vector{Float64}     # P-C coordinates
    FnT::AbstractVector     # Boys function values
    mtp::AbstractVector     # (-2p)^n powers
end 


function base(t, u, v, n, pair::HermitePair)
    #  n             n
    # R      = (-2 p)  F (T)
    #  0,0,0            n

    Rtuv = pair.FnT[n+1] * pair.mtp[n+1]

    Rtuv
end


function vertical_t(t, u, v, n, pair::HermitePair)
    #  n              n+1            n+1
    # R      = (t-1) R        + X   R
    #  t,u,v          t-2,u,v    PC  t-1,u,v

    Rtuv = hermite(t-1, u, v, n+1, pair) * pair.PC[1]
    if t > 1
        Rtuv += hermite(t-2, u, v, n+1, pair) * (t-1)
    end

    Rtuv
end


function vertical_u(t, u, v, n, pair::HermitePair)
    #  n              n+1            n+1
    # R      = (u-1) R        + Y   R
    #  t,u,v          t,u-2,v    PC  t,u-1,v

    Rtuv = hermite(t, u-1, v, n+1, pair) * pair.PC[2]
    if u > 1
        Rtuv += hermite(t, u-2, v, n+1, pair) * (u-1)
    end

    Rtuv
end


function vertical_v(t, u, v, n, pair::HermitePair)
    #  n              n+1            n+1
    # R      = (v-1) R        + Z   R
    #  t,u,v          t,u,v-2    PC  t,u,v-1
    
    Rtuv = hermite(t, u, v-1, n+1, pair) * pair.PC[3]
    if v > 1
        Rtuv += hermite(t, u, v-2, n+1, pair) * (v-1)
    end

    Rtuv
end


function get_hermite(t, u, v, n, pair::HermitePair)
    # Default path for computing terms.

    if (t < 0) || (u < 0) || (v < 0)
        Rtuv = 0.0
    elseif t == u == v == 0
        Rtuv = base(t, u, v, n, pair)
    elseif t == u == 0
        Rtuv = vertical_v(t, u, v, n, pair)
    elseif t == 0
        Rtuv = vertical_u(t, u, v, n, pair)
    else
        Rtuv = vertical_t(t, u, v, n, pair)
    end

    Rtuv
end