# Line algorithms

"""For convinience"""
dp_line_naive(y::Vector{ℝ}, λ::ℝ) = dp_line_naive(y, i->λ, i->1.0)

"""For convinience..."""
dp_line(y::Vector{ℝ}, λ::ℝ) = dp_line(y, i->λ, i->1.0)


"""Naive implementation (based on PWLs) for path graphs, i.e. 1d fused LASSO"""
function dp_line_naive(y, λ, μ)
    n = length(y)
    lb, ub = fill(-∞, n), fill(+∞, n)
    q(i) = PWL(0.0, -μ(i)*y[i]; slope=μ(i))
    ∂f = q(1)
    for i = 1:(n-1)
        lb[i] = find_x(∂f, -λ(i))
        ub[i] = find_x(∂f, +λ(i))
        ∂f = q(i+1) + clip(∂f, -λ(i), +λ(i))
    end
    xn = find_x(∂f, 0)
    @debug "\nlb = $lb\nub = $ub"
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
function dp_line(y, λ, µ)
    n = length(y)
    lb, ub = fill(∞, n), fill(-∞, n)
    σ(i) = if i > 1; λ(i-1) else 0.0 end
    line(i, r) = LineSegment(µ(i), -µ(i)*y[i] + r)
    pq = DeQue{Event}()
    for i = 1:(n-1)
        lb[i] = clip_front(pq, line(i, -σ(i)), -λ(i))
        ub[i] = clip_back(pq,  line(i, +σ(i)), +λ(i))
    end
    @debug "lb = $lb"
    @debug "ub = $ub"
    @debug "pq = $pq"
    xn = clip_front(pq, line(n, -σ(n)), 0.0)
    return dp_line_backtrace(xn, lb, ub)
end
