"""Line algorithms"""

function dp_line_naive(y::Vector{Float64}, lambda::Float64)
    n = length(y)
    lb, ub, x = fill(Inf, n), fill(-Inf, n), zeros(n)
    q(i) = PWL(0.0, -y[1]; slope=1.0)
    df = q(1)
    for i = 2:n
        lb[i-1] = find_x(df, -lambda)
        ub[i-1] = find_x(df, +lambda)
        df = q(i) + clip(df, lb[i-1], ub[i-1])
    end

    x[n] = find_x(df, 0)
    for i = n:-1:2
        x[i-1] = clip(x[i], lb[i-1], ub[i-1])
    end

    return x
end

        
