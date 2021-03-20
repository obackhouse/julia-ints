# Compute the expansion coefficients

# KAB::Float64    # exponential term
# PA::Float64     # P-A coordinate
# PB::Float64     # P-B coordinate
# p::Float64      # α+β
# q::Float64      # 1/(2(α+β))


function base(i, j, t, KAB, PA, PB, p, q)
    #  0,0
    # E    = K
    #  0      AB
    
    Eab = KAB
    
    Eab
end


function oob(i, j, t, KAB, PA, PB, p, q)
    #  i,j
    # E    = 0.0  if t<0 or i<0 or j<0 or t>(i+j)
    #  t
    
    Eab = zero(KAB)

    Eab
end


function mcmurchie_davidson_i(i, j, t, KAB, PA, PB, p, q)
    #  i,j     1    i-1,j        i-1,j          i-1,j
    # E    = ----- E      + X   E      + (t+1) E
    #  t      2*p   t-1      PA  t              t+1

    Eab = expansion(i-1, j, t-1, KAB, PA, PB, p, q) .* q .+
          expansion(i-1, j, t,   KAB, PA, PB, p, q) .* PA .+
          expansion(i-1, j, t+1, KAB, PA, PB, p, q) * (t+1)

    Eab
end


function mcmurchie_davidson_j(i, j, t, KAB, PA, PB, p, q)
    #  i,j     1    i,j-1        i,j-1          i,j-1
    # E    = ----- E      + X   E      + (t+1) E
    #  t      2*p   t-1      PA  t              t+1

    Eab = expansion(i, j-1, t-1, KAB, PA, PB, p, q) .* q .+
          expansion(i, j-1, t,   KAB, PA, PB, p, q) .* PB .+
          expansion(i, j-1, t+1, KAB, PA, PB, p, q) * (t+1)

    Eab
end


function two_term(i, j, t, KAB, PA, PB, p, q)
    #  i,j     1   /    i-1,j      i,j-1 \
    # E    = ----- | i E      + j E      |   for t > 0
    #  t      2pt  \    t-1        t-1   /
    
    Eab = (expansion(i-1, j,   t-1, KAB, PA, PB, p, q) * i .+
           expansion(i,   j-1, t-1, KAB, PA, PB, p, q) * j) .* (q / t)

    Eab
end


function order_dump_i(i, j, t, KAB, PA, PB, p, q)
    #  i,0   /   1   \ t / i \  i-t,0
    # E    = | ----- |   |   | E
    #  t     \  2*p  /   \ t /  0

    Eab = expansion(i-t, 0, 0, KAB, PA, PB, p, q) .* q.^t * binomial(i, t)

    Eab
end


function order_dump_j(i, j, t, KAB, PA, PB, p, q)
    #  0,j   /   1   \ t / j \  0,j-t
    # E    = | ----- |   |   | E
    #  t     \  2*p  /   \ t /  0

    Eab = expansion(0, j-t, 0, KAB, PA, PB, p, q) .* q.^t * binomial(j, t)

    Eab
end


function obara_saika_i(i, j, t, KAB, PA, PB, p, q)
    #  i,j        i-1,j     1   /        i-2,j      i-1,j-1    i-1,j \
    # E    = X   E      + ----- | (i-1) E      + j E        + E      |
    #  t      PA  t        2*p  \        t          t          t-1   /

    Eab = expansion(i-1, j,   t,   KAB, PA, PB, p, q) .* PA .+ q .* (
          expansion(i-2, j,   t,   KAB, PA, PB, p, q) * (i-1) .+
          expansion(i-1, j-1, t,   KAB, PA, PB, p, q) * j .+
          expansion(i-1, j,   t-1, KAB, PA, PB, p, q)
    )

    Eab
end


function obara_saika_j(i, j, t, KAB, PA, PB, p, q)
    #  i,j        i,j-1     1   /    i-1,j-1          i,j-2    i,j-1 \
    # E    = X   E      + ----- | i E        + (j-1) E      + E      |
    #  t      PB  t        2*p  \    t                t        t-1   /
    
    Eab = expansion(i,   j-1, t,   KAB, PA, PB, p, q) .* PB .+ q .* (
          expansion(i-1, j-1, t,   KAB, PA, PB, p, q) * i .+
          expansion(i,   j-2, t,   KAB, PA, PB, p, q) * (j-1) .+
          expansion(i,   j-1, t-1, KAB, PA, PB, p, q)
    )

    Eab
end


function expansion(i, j, t, KAB, PA, PB, p, q)
    # Default path for computing terms - doesn't use all of the
    # recurrences and better expansions probably exist.

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


function expansion_mmd(i, j, t, KAB, PA, PB, p, q)
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
