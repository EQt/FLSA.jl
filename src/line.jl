"""Line algorithms"""

function dp_line_naive(y::Vector{Float64}, lambda::Float64)
    n = length(y)
    lb, ub, x = fill(Inf, n), fill(-Inf, n), zeros(n)
    df = PWL(0.0, -y[1]; slope=1.0)
    for i = 2:n
        lb[i-1] = find_x(df, -lambda)
        ub[i-1] = find_x(df, +lambda)
        df = PWL(0.0, -y[i]; slope=1.0) + clip(df, lb[i-1], ub[i-1])
    end

    x[n] = find_x(df, 0)
    for i = n:-1:2
        x[i-1] = clip(x[i], lb[i-1], ub[i-1])
    end

    return x
end

        
