"""Line algorithms"""

function dp_line_naive(y::Vector{Float64}, λ::Float64)
    n = length(y)
    lb, ub = fill(Inf, n), fill(-Inf, n)

    q(i) = PWL(0.0, -y[i]; slope=1.0)
    df = q(1)
    for i = 2:n
        lb[i-1] = find_x(df, -λ)
        ub[i-1] = find_x(df, +λ)
        df = q(i) + clip(df, lb[i-1], ub[i-1])
    end

    xn = find_x(df, 0)
    return dp_line_backtrace(xn, lb, ub)
end


function dp_line_backtrace(xn, lb, ub)
    n = length(lb)
    x = zeros(n)
    x[n] = xn
    for i = n:-1:2
        x[i-1] = clip(x[i], lb[i-1], ub[i-1])
    end
    return x
end


function dp_line(y::Vector{Float64}, λ::Float64)
    n = length(y)
    lb, ub = fill(Inf, n), fill(-Inf, n)
    pq = Deque{Event}()
    o1, o2 = 0.0, 0.0 # TODO
    push_front!(pq, Event(y[1]-λ, o1, +1.0))
    push_back!(pq,  Event(y[1]+λ, o2, -1.0))
    find_min(pq, x) = x # TODO
    find_max(pq, x) = x # TODO
    for i = 2:n
        lb[i-1] = find_min(pq, -λ)
        ub[i-1] = find_max(pq, +λ)
    end

    xn = find_min(pq, 0)
    return dp_line_backtrace(xn, lb, ub)
end
