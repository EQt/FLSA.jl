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
    child(i) = if i > 1;  -λ(i-1) else 0.0 end
    line_f(i) = LineSegment(µ(i), -µ(i)*y[i] - child(i))
    line_b(i) = LineSegment(µ(i), -µ(i)*y[i] + child(i))
    pq = DeQue{Event}()
    for i = 1:(n-1)
        lb[i] = clip_front(pq, init_front(i), -λ(i))
        ub[i] = clip_back(pq,  init_back(i),  +λ(i))
    end
    @debug "lb = $lb"
    @debug "ub = $ub"
    @debug "pq = $pq"
    xn = clip_front(pq, init_front(n), 0.0)
    return dp_line_backtrace(xn, lb, ub)
end


"""Clip PWL represented by `pq` from negative until `t`, starting with `l`"""
function clip_front{Q}(pq::Q, l::LineSegment, t::ℝ)
    x = find_x(t, l)
    @debug "Starting with l=$l, x=$x, min_x = $(min_x(pq))"
    while x > min_x(pq)
        e = pop_front!(pq)
        @debug "pop_front  pq = $pq, e=$e"
        l.offset += e.offset
        l.slope  += e.slope
        x = find_x(t, l)
    end
    push_front!(pq, Event(x, l.slope, l.offset - t))
    @debug "push_front pq = $pq"
    return x
end


"""Same as `clip_front` just everything reversed"""
function clip_back{Q}(pq::Q, l::LineSegment, t::ℝ)
    x = find_x(t, l)
    while max_x(pq) > x
        e = pop_back!(pq)
        l.offset -= e.offset
        l.slope  -= e.slope
        x = find_x(t, l)
    end
    push_back!(pq, Event(x, -l.slope, -l.offset + t))
    @debug "push_back pq = $pq"
    return x
end


"""Extract minimal position `x` of queue `pq` or `+∞` if none exists"""
@inline min_x{Q}(pq::Q) = try front(pq).x catch +∞ end

"""Extract maximal position `x` of queue `pq` or `-∞` if none exists"""
@inline max_x{Q}(pq::Q) = try back(pq).x  catch -∞ end


"""For convinience..."""
dp_line(y::Vector{Float64}, λ::Float64) = dp_line(y, i->λ, i->1.0)
