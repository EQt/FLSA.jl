# Line algorithms

"""Naive implementation (based on PWLs) for path graphs, i.e. 1d fused LASSO"""
function dp_line_naive(y, λ, μ)
    n = length(y)
    lb, ub = fill(∞, n), fill(-∞, n)
    q(i) = PWL(0.0, -μ(i)*y[i]; slope=μ(i))
    ∂f = q(1)
    for i = 1:(n-1)
        lb[i] = find_x(∂f, -λ(i))
        ub[i] = find_x(∂f, +λ(i))
        ∂f = q(i+1) + clip(∂f, -λ(i), +λ(i))
    end
    xn = find_x(∂f, 0)
    return dp_line_backtrace(xn, lb, ub)
end


"""For convinience"""
dp_line_naive(y::Vector{Float64}, λ::Float64) = dp_line_naive(y, i->λ, i->1.0)


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
    init_up(i)   = LineSegment(µ(i), -λ(i) - µ(i)*y[i])
    init_down(i) = LineSegment(µ(i), +λ(i) - µ(i)*y[i])
    pq = DeQue{Event}()
    for i = 1:(n-1)
        lb[i] = clip_up(pq,   init_up(i),   -λ(i))
        ub[i] = clip_down(pq, init_down(i), +λ(i))
    end
    @debug "lb = $lb"
    @debug "ub = $ub"
    xn = clip_up(pq, init_up(n), 0.0)
    return dp_line_backtrace(xn, lb, ub)
end


"""Clip PWL represented by `pq` from negative until `t`, starting with `l`"""
function clip_up{Q}(pq::Q, l::LineSegment, t::ℝ)
    @debug "Starting with l=$l"
    x = find_x(t, l)
    while min_x(pq) ≤ x
        e = pop_front!(pq)
        l.offset += e.offset
        l.slope  += e.slope
        x = find_x(t, l)
    end
    push_front!(pq, Event(x, l))
    return x
end


"""Same as `clip_up` just everything reversed"""
function clip_down{Q}(pq::Q, l::LineSegment, t::ℝ)
    x = find_x(t, l)
    while x ≤ max_x(pq)
        e = pop_back!(pq)
        l.offset -= e.offset
        l.slope  -= e.slope
        x = find_x(t, l)
    end
    push_back!(pq, Event(x, l))
    return x
end


"""Extract position of minimal `x` or `∞` if none exists"""
@inline min_x{Q}(pq::Q) = try front(pq).x catch +∞ end
@inline max_x{Q}(pq::Q) = try back(pq).x  catch -∞ end


"""For convinience..."""
dp_line(y::Vector{Float64}, λ::Float64) = dp_line(y, i->λ, i->1.0)
