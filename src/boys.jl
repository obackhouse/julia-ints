# Evaluate the Boys function

include("boys_lookup.jl")
include("utils.jl")


function boys_long(n, T)
    # Compute the Boys function Fn(T) with a long-range approximation

    fac = factorial_2(2*n-1)
    phase = 2 ^ (-1-n)
    Tpow = T ^ (-0.5-n)
    FnT = fac * √π * phase * Tpow

    FnT
end


function boys(n, T)
    # Compute the Boys function Fn(T)

    if T == 0
        FnT = 1.0 / (2 * n + 1)
    elseif T < BOYS_LONG_RANGE_TMAX
        idx = round(Int64, T * BOYS_RESOLUTION) + 1
        
        dT = (idx - T * BOYS_RESOLUTION) * (1.0 / BOYS_RESOLUTION)
        fac = 1.0
        FnT = 0.0

        for n = 1:(BOYS_INTERP_NUM-1)
            fac *= (dT / n)
            FnT += fac * boys_fn_data[n+i+1][idx]
        end
    else
        FnT = boys_long(n, T)
    end

    FnT
end


function boys_array!(nmax, T, out::AbstractVector)
    # Populate an array with the values of the Boys function Fn(T) up
    # to order nmax
    
    @assert nmax >= 0
    @assert T >= 0
    @assert nmax <= BOYS_MAX_DATA

    if T == 0
        for n = 0:nmax
            out[n] = 1.0 / (2 * n + 1)
        end
    elseif T < BOYS_LONG_RANGE_TMAX
        idx = round(Int64, T * BOYS_RESOLUTION) + 1
        out[nmax+1] = boys_fn_data[nmax+1][idx]

        dT = (idx - T * BOYS_RESOLUTION) * (1.0 / BOYS_RESOLUTION)
        fac = 1.0

        for n = 1:(BOYS_INTERP_NUM-1)
            fac *= (dT / n)
            out[nmax] += fac * boys_fn_data[nmax+n+1][idx]
        end
    else
        out[nmax+1] = boys_long(nmax, T)
    end

    if nmax > 0
        expT = exp(-1.0 * T)
        
        for n = nmax:1
            out[n] = (2 * T * out[n] + expT) / (2 * n - 1)
        end
    end

    out
end
