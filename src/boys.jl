# Evaluate the Boys function

include("boys_lookup.jl")
include("utils.jl")


function boys_long(n::Int64, T::Float64)
    # Compute the Boys function Fn(T) with a long-range approximation

    fac = factorial_2(2*n-1)
    phase = 2.0^(-1-n)
    Tpow = T^(-0.5-n)
    FnT = fac * √π * phase * Tpow

    FnT
end


function boys(n::Int64, T::Float64)
    # Compute the Boys function Fn(T)

    if T < BOYS_TZERO
        FnT = 1.0 / (2 * n + 1)
    elseif T < BOYS_LONG_RANGE_TMAX
        idx = round(Int64, T * BOYS_RESOLUTION)
        
        dT = (idx - T * BOYS_RESOLUTION) * (1.0 / BOYS_RESOLUTION)
        fac = 1.0
        FnT = 0.0

        for i = 1:(BOYS_INTERP_NUM-1)
            fac *= (dT / i)
            FnT += fac * boys_fn_data[n+i+1][idx+1]
        end
    else
        FnT = boys_long(n, T)
    end

    FnT
end


function boys_array!(nmax::Int64, T::Float64, out::AbstractVector{Float64})
    # Populate an array with the values of the Boys function Fn(T) up
    # to order nmax
    
    @assert nmax >= 0
    @assert T >= 0
    @assert nmax <= BOYS_MAX_DATA

    if T < BOYS_TZERO
        for n = 0:nmax
            out[n+1] = 1.0 / (2 * n + 1)
        end
    else
        if T < BOYS_LONG_RANGE_TMAX
            idx = round(Int64, T * BOYS_RESOLUTION)
            out[nmax+1] = boys_fn_data[nmax+1][idx+1]

            dT = (idx - T * BOYS_RESOLUTION) / BOYS_RESOLUTION
            fac = 1.0

            for n = 1:(BOYS_INTERP_NUM-1)
                fac *= (dT / n)
                out[nmax+1] += fac * boys_fn_data[nmax+n+1][idx+1]
            end
        else
            out[nmax+1] = boys_long(nmax, T)
        end

        if nmax > 0
            expT = exp(-1.0 * T)
            
            for n = nmax:-1:1
                out[n] = (2 * T * out[n+1] + expT) / (2 * n - 1)
            end
        end
    end

    out
end
