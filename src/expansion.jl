# Compute the expansion coefficients
#

struct ExpansionPair
    Kab::Float64    # exponential term
    Pa::Float64     # P-A coordinate
    Pb::Float64     # P-B coordinate
    α::Float64      # α exponent
    β::Float64      # β exponent
    p::Float64      # α+β
    q::Float64      # 1/(2(α+β))
end


function base(i, j, t, pair::ExpansionPair)
    #  0,0
    # E    = K
    #  0      AB
    
    Eab = pair.Kab
    
    Eab
end


function mcmurchie_davidson_i(i, j, t, pair::ExpansionPair)
    #  i,j     1    i-1,j        i-1,j          i-1,j
    # E    = ----- E      + X   E      + (t+1) E
    #  t      2*p   t-1      PA  t              t+1
    
    Eab  = expansion(i-1, j, t-1, pair) * pair.q
    Eab += expansion(i-1, j, t,   pair) * pair.β
    Eab += expansion(i, j-1, t+1, pair) * (t+1)

    Eab
end


function mcmurchie_davidson_j(i, j, t, pair::ExpansionPair)
    #  i,j     1    i,j-1        i,j-1          i,j-1
    # E    = ----- E      + X   E      + (t+1) E
    #  t      2*p   t-1      PA  t              t+1

    Eab  = expansion(i, j-1, t-1, pair) * pair.q
    Eab += expansion(i, j-1, t,   pair) * pair.α
    Eab += expansion(i, j-1, t+1, pair) * (t+1)

    Eab
end


function two_term(i, j, t, pair::ExpansionPair)
    #  i,j     1   /    i-1,j      i,j-1 \
    # E    = ----- | i E      + j E      |   for t > 0
    #  t      2pt  \    t-1        t-1   /
    
    Eab = 0.0

    if i > 0
        Eab += expansion(i-1, j, t-1, pair) * i
    end
    if j > 0
        Eab += expansion(i, j-1, t-1, pair) * j
    end
    Eab *= (pair.q / t)

    Eab
end


function order_dump_i(i, j, t, pair::ExpansionPair)
    #  i,0   /   1   \ t / i \  i-t,0
    # E    = | ----- |   |   | E
    #  t     \  2*p  /   \ t /  0

    Eab  = expansion(i-t, 0, 0, pair)
    Eab *= binomial(i, t)
    Eab *= pair.q^t

    Eab
end


function order_dump_j(i, j, t, pair::ExpansionPair)
    #  0,j   /   1   \ t / i \  0,j-t
    # E    = | ----- |   |   | E
    #  t     \  2*p  /   \ j /  0

    Eab  = expansion(0, j-t, 0, pair)
    Eab *= binomial(j, t)
    Eab *= pair.q^t

    Eab
end


function obara_saika_i(i, j, t, pair::ExpansionPair)
    #  i,j        i-1,j     1   /        i-2,j      i-1,j-1    i-1,j \
    # E    = X   E      + ----- | (i-1) E      + j E        + E      |
    #  t      PA  t        2*p  \        t          t          t-1   /
    
    Eab = 0.0

    if i > 1
        Eab += expansion(i-2, j, t, pair) * (i-1)
    end
    if j > 0
        Eab += expansion(i-1, j-1, t, pair) * j
    end
    Eab += expansion(i-1, j, t-1, pair)
    Eab *= pair.q
    Eab += expansion(i-1, j, t, pair) * pair.α

    Eab
end


function obara_saika_j(i, j, t, pair::ExpansionPair)
    #  i,j        i,j-1     1   /    i-1,j-1          i,j-2    i,j-1 \
    # E    = X   E      + ----- | i E        + (j-1) E      + E      |
    #  t      PB  t        2*p  \    t                t        t-1   /

    if i > 0
        Eab += expansion(i-1, j-1, t, pair) * i
    end
    if j > 1
        Eab += expansion(i, j-2, t, pair) * (j-1)
    end
    Eab += expansion(i, j-1, t-1, pair)
    Eab *= pair.q
    Eab += expansion(i, j-1, t, pair) * pair.β

    Eab
end


function expansion(i, j, t, pair::ExpansionPair)
    # Default path for computing terms - doesn't use all of the
    # recurrences and better expansions probably exist.

    if (t < 0) || (i < 0) || (j < 0) || (t > (i+j))
        Eab = 0.0
    elseif i == j == t == 0
        Eab = base(i, j, t, pair)
    elseif j == 0 && t > 0
        Eab = order_dump_i(i, j, t, pair)
    elseif i == 0 && t > 0
        Eab = order_dump_j(i, j, t, pair)
    elseif j == 0
        Eab = obara_saika_i(i, j, t, pair)
    else
        Eab = obara_saika_j(i, j, t, pair)
    end

    Eab
end
