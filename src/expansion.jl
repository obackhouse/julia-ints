# Compute the expansion coefficients

# KAB::Float64    # exponential term
# PA::Float64     # P-A coordinate
# PB::Float64     # P-B coordinate
# p::Float64      # α+β
# q::Float64      # 1/(2(α+β))


function base(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  0,0
    # E    = K
    #  0      AB
    
    KAB
end


function oob(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,j
    # E    = 0.0  if t<0 or i<0 or j<0 or t>(i+j)
    #  t
    
    zero(KAB)
end


function mcmurchie_davidson_i(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,j     1    i-1,j        i-1,j          i-1,j
    # E    = ----- E      + X   E      + (t+1) E
    #  t      2*p   t-1      PA  t              t+1

    Eab = expansion(i-1, j, t-1, KAB, PA, PB, p, q) * q +
          expansion(i-1, j, t,   KAB, PA, PB, p, q) * PA +
          expansion(i-1, j, t+1, KAB, PA, PB, p, q) * (t+1)

    Eab
end


function mcmurchie_davidson_j(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,j     1    i,j-1        i,j-1          i,j-1
    # E    = ----- E      + X   E      + (t+1) E
    #  t      2*p   t-1      PA  t              t+1

    Eab = expansion(i, j-1, t-1, KAB, PA, PB, p, q) * q +
          expansion(i, j-1, t,   KAB, PA, PB, p, q) * PB +
          expansion(i, j-1, t+1, KAB, PA, PB, p, q) * (t+1)

    Eab
end


function two_term(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,j     1   /    i-1,j      i,j-1 \
    # E    = ----- | i E      + j E      |   for t > 0
    #  t      2pt  \    t-1        t-1   /
    
    Eab = (expansion(i-1, j,   t-1, KAB, PA, PB, p, q) * i +
           expansion(i,   j-1, t-1, KAB, PA, PB, p, q) * j) * (q / t)

    Eab
end


function order_dump_i(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,0   /   1   \ t / i \  i-t,0
    # E    = | ----- |   |   | E
    #  t     \  2*p  /   \ t /  0

    Eab = expansion(i-t, 0, 0, KAB, PA, PB, p, q) * q^t * binomial(i, t)

    Eab
end


function order_dump_j(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  0,j   /   1   \ t / j \  0,j-t
    # E    = | ----- |   |   | E
    #  t     \  2*p  /   \ t /  0

    Eab = expansion(0, j-t, 0, KAB, PA, PB, p, q) * q^t * binomial(j, t)

    Eab
end


function obara_saika_i(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,j        i-1,j     1   /        i-2,j      i-1,j-1    i-1,j \
    # E    = X   E      + ----- | (i-1) E      + j E        + E      |
    #  t      PA  t        2*p  \        t          t          t-1   /

    Eab = expansion(i-1, j,   t,   KAB, PA, PB, p, q) * PA + q * (
          expansion(i-2, j,   t,   KAB, PA, PB, p, q) * (i-1) +
          expansion(i-1, j-1, t,   KAB, PA, PB, p, q) * j +
          expansion(i-1, j,   t-1, KAB, PA, PB, p, q)
    )

    Eab
end


function obara_saika_j(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    #  i,j        i,j-1     1   /    i-1,j-1          i,j-2    i,j-1 \
    # E    = X   E      + ----- | i E        + (j-1) E      + E      |
    #  t      PB  t        2*p  \    t                t        t-1   /
    
    Eab = expansion(i,   j-1, t,   KAB, PA, PB, p, q) * PB + q * (
          expansion(i-1, j-1, t,   KAB, PA, PB, p, q) * i +
          expansion(i,   j-2, t,   KAB, PA, PB, p, q) * (j-1) +
          expansion(i,   j-1, t-1, KAB, PA, PB, p, q)
    )

    Eab
end


function expansion_os(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    # Obara-Saika path with some high-order optimisations

    if (t < 0) || (i < 0) || (j < 0) || (t > (i+j))
        Eab = oob(i, j, t, KAB, PA, PB, p, q)
    elseif i == j == t == 0
        Eab = base(i, j, t, KAB, PA, PB, p, q)
    elseif j == 0 && t > 0
        Eab = order_dump_i(i, j, t, KAB, PA, PB, p, q)
    elseif i == 0 && t > 0
        Eab = order_dump_j(i, j, t, KAB, PA, PB, p, q)
    elseif j == 0
        Eab = obara_saika_i(i, j, t, KAB, PA, PB, p, q)
    else
        Eab = obara_saika_j(i, j, t, KAB, PA, PB, p, q)
    end

    Eab
end


function expansion_mmd(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
)
    # Standard McMurchie-Davidson path

    if (t < 0) || (t > (i+j)) || (i < 0) || (j < 0)
        Eab = oob(i, j, t, KAB, PA, PB, p, q)
    elseif i == j == t == 0
        Eab = base(i, j, t, KAB, PA, PB, p, q)
    elseif j == 0
        Eab = mcmurchie_davidson_i(i, j, t, KAB, PA, PB, p, q)
    else
        Eab = mcmurchie_davidson_j(i, j, t, KAB, PA, PB, p, q)
    end

    Eab
end


function expansion(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
        cache::Nothing=nothing,
)
    # Driver function for expansion coefficients

    expansion_os(i, j, t, KAB, PA, PB, p, q)
end


function expansion(
        i::Int64, 
        j::Int64, 
        t::Int64, 
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
        cache::Dict{NTuple{3, Int64}, Float64},
)
    # Driver function for expansion coefficients with caching

    key = (i, j, t)
    if haskey(cache, key)
        return cache[key]
    end

    Eab = expansion_os(i, j, t, KAB, PA, PB, p, q)

    push!(cache, key => Eab)

    Eab
end


function populate_expansion(
        i::Int64,
        j::Int64,
        tmax::Int64,
        KAB::Float64, 
        PA::Float64, 
        PB::Float64, 
        p::Float64, 
        q::Float64, 
        out::AbstractArray{Float64},
)
    # Populate array with expansion coefficients for t = 0 → tmax

    if (i+j+tmax) > 5  #TODO: I assume this is good for high ang mom?
        cache = Dict{NTuple{3, Int64}, Float64}()
    else
        cache = nothing
    end

    #TODO more efficient path
    for t = 0:tmax
        out[t+1] = expansion(i, j, t, KAB, PA, PB, p, q, cache)
    end

    out
end
