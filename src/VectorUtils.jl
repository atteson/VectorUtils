module VectorUtils

export brent

function brent( f, x0::AbstractVector{T}, x1::AbstractVector{T};
                EPS = eps(T),
                xtol = 1e-7,
                ytol = 2*EPS,
                maxiter::Integer = 50,
                ) where {T}
    y0 = f(x0)
    y1 = f(x1)

    ay0lay1 = abs.(y0) .< abs.(y1)
    
    tmp = x0
    x0 = ifelse.( ay0lay1, x1, x0 )
    x1 = ifelse.( ay0lay1, tmp, x1 )
    
    tmp = y0
    y0 = ifelse.( ay0lay1, y1, y0 )
    y1 = ifelse.( ay0lay1, tmp, y1 )
    
    x2 = copy(x0)
    y2 = copy(y0)
    x3 = copy(x2)

    bisection = .!isnan.(x0)
    for _ in 1:maxiter
        if all(abs.(x1 - x0) .< xtol)
            return (x1,y1)
        end

        x = ifelse.(
            (abs.(y0-y2) .> ytol) .& (abs.(y1-y2) .> ytol),
            
            x0 .* y1 .* y2 ./((y0 - y1) .* (y0 - y2))
            + x1 .* y0 .* y2 ./((y1 - y0) .* (y1 - y2))
            + x2 .* y0 .* y1 ./((y2 - y0) .* (y2 - y1)),
            
            x1 - y1 .* (x1 - x0)./(y1 - y0)
        )

        delta = 2 * EPS * abs.(x1)
        min1 = abs.(x - x1)
        min2 = abs.(x1 - x2)
        min3 = abs.(x2 - x3)

        bisection = (
            ((x .< (3 .* x0 + x1)./4) .| (x .> x1)) .|
            ((bisection .!=0 ) .& (min1 .>= min2./2)) .|
            (.!bisection .& (min1 .>= min3/2)) .|
            (bisection .& (min2 .< delta)) .|
            (.!bisection .& (min3 .< delta))
        )

        x = ifelse.( bisection, (x0 + x1)./2, x )

        y = f(x)
        if all(abs.(y) .< ytol)
            return (x,y)
        end
        x3 = x2
        x2 = x1
        signs = sign.(y0) .!= sign.(y)
        x1 = ifelse.( signs, x, x1 )
        y1 = ifelse.( signs, y, y1 )
        x0 = ifelse.( .!signs, x, x0 )
        y1 = ifelse.( .!signs, y, y0 )
    end
    error( "Maximum iterations exceeded!" )
end

end # module VectorUtils
