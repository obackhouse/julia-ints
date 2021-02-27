# Utility functions


function factorial_2(n)
    # Compute the double factorial n!!

    if (n == 1 || n == 0 || n == -1)
        fac = 1
    else
        fac = n * factorial_2(n-2)
    end

    fac
end
