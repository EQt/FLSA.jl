# Line algorithms

"""Naive implementation (based on PWLs) for path graphs, i.e. 1d fused LASSO"""
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

"""FLSA on a line, computed by Johnson's fast *dynamic programming* algorithm"""
function dp_line(y::Vector{Float64}, λ::Float64)
    # TODO: compare to file ../experimental/treepwl_error1.jl
    n = length(y)
    lb, ub = fill(Inf, n), fill(-Inf, n)
    pq = Deque{Event}()
    o1, o2 = 0.0, 0.0 # TODO
    push_front!(pq, Event(y[1]-λ, o1, +1.0))
    push_back!(pq,  Event(y[1]+λ, o2, -1.0))
    find_min(pq, x) = x # TODO
    find_max(pq, x) = x # TODO
    for i = 2:n
        lb[i-1] = min_event(pq, -λ)
        ub[i-1] = max_event(pq, +λ)
    end

    xn = find_min(pq, 0)
    return dp_line_backtrace(xn, lb, ub)
end


"""Position where the line with given `offset` and `slope` will have value `c`"""
@inline find_x(offset, slope, c) = (c - offset)/slope

"""Extract position of minimal `x` or `∞` if none exists"""
@inline min_x{Q}(pq::Q) = try front(pq).x catch throw(∞) end

"""Compute a new lower bound event for node `v`"""
function min_event(pq, v::Int, c::Float64, y)
    local slope = 1.0
    local offset = y[v] - c
    local x = find_x(offset, slope, c)
    xk = min_x(pq)
    while x > xk
        e = pop_front!(pq)
        offset += e.offset
        slope  += e.slope
        xk = min_x(pq)
        x = forecast(offset, slope, c)
    end
    push_front!(pq, Event(x, offset, slope))
    return x
end
