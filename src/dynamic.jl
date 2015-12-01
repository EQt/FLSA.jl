import Base.string


"""
After computing the bounds [lb, ub] for each node, compute optimal solution `x`
by clipping each edge, i.e. backtracing from root value `xr` down to children
"""
function dp_tree_backtrace(xr::ℝ, t::Tree, lb, ub)
    @debug "backtrace: xr=$xr"
    x = zeros(num_nodes(t))
    x[t.root] = xr
    for v in preorder(t)
        x[v] = clip(x[t.parent[v]], lb[v], ub[v])
        @debug "x[$v] = $(x[v]), x[parent] = $(x[t.parent[v]])"
    end
    return x
end


"""Compute dual solution α on a tree, such that y = D' α"""
function dual_tree(y::Vector{ℝ}, tree)
    α = zeros(length(tree.edges))
    for c in postorder(tree)
        v = tree.parent[c]
        e = edge_index(tree, (c,v))
        α[e]  = sign(c-v) * y[c]
        y[v] += sign(c-v) * α[e]
    end
    return α
end

"""For convinience..."""
dual_tree(y::Vector{Float64}, x::Vector{Float64}, tree) = dual_tree(x-y, tree)


"""Represent a line segment"""
type LineSegment
    slope::ℝ
    offset::ℝ
    function LineSegment(s, o)
        @assert abs(s) > 1e-15
        new(s,o)
    end
end


"""Record what is happening, when a knot of the PWL is hit"""
immutable Event
    x::Float64      # position
    slope::Float64  # delta slope
    offset::Float64 # delta offset
    function Event(x, s, o)
        @assert isfinite(x)
        @assert abs(s) > 1e-16
        new(x, s, o)
    end
end

@inline event_time(e::Event) = e.x

include("queues.jl")

"""Print them more readable"""
string(e::Event) = "$(e.x) @ $(e.slope)x + $(e.slope))"


"""Find x, such that t = slope*x + offset"""
@inline find_x(t::ℝ, slope::ℝ, offset::ℝ) = (t - offset)/slope

@inline function find_x(t::ℝ, ls::LineSegment)
    x = find_x(t, ls.slope, ls.offset)
    @assert isfinite(x) "t=$t, ls=$ls"
    x
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
@inline min_x{Q}(pq::Q) = if !isempty(pq); front(pq).x else +∞ end

"""Extract maximal position `x` of queue `pq` or `-∞` if none exists"""
@inline max_x{Q}(pq::Q) = if !isempty(pq); back(pq).x  else -∞ end


"""For convinience..."""
dp_tree_naive(y::Vector{ℝ}, λ::ℝ, t::Tree) = dp_tree_naive(y, i->λ, i->1.0, t)

"""Compute x=FLSA(y, λ) on a a (sub)tree t, naive PWL implementation"""
function dp_tree_naive(y::Vector{ℝ}, λ, µ, t::Tree)
    n = length(y)
    ∂f = [PWL(0.0, -µ(i)*y[i]; slope=µ(i)) for i=1:n]
    lb, ub = fill(-∞, n), fill(∞, n)
    for i in postorder(t)
        @debug "processing node $i: (n=$n)"
        lb[i] = find_x(∂f[i], -λ(i))
        ub[i] = find_x(∂f[i], +λ(i))
        ∂f[t.parent[i]] += clip(∂f[i], -λ(i), +λ(i))
    end
    @debug @val lb
    @debug @val ub
    @debug @val ∂f[t.root]
    xr = find_x(∂f[t.root], 0.0)
    dp_tree_backtrace(xr, t, lb, ub)
end


"""For convinience..."""
dp_tree(y::Vector{ℝ}, λ::ℝ, t::Tree) = dp_tree(y, i->λ, i->1.0, t)

"""FLSA on a line, computed by Johnson's fast *dynamic programming* algorithm"""
function dp_tree(y, λ, µ, t::Tree)
    n, r = length(y), t.root
    lb, ub = fill(-∞, n), fill(∞, n)
    σ(i) = sum(ℝ[λ(c) for c in t.children[i]])
    line(i, r) = LineSegment(µ(i), -µ(i)*y[i] + r)
    pq = [EventQueue() for i = 1:n]
    for i in postorder(t)
        lb[i] = clip_front(pq[i], line(i, -σ(i)), -λ(i))
        ub[i] = clip_back( pq[i], line(i, +σ(i)), +λ(i))
        pq[t.parent[i]] += pq[i]
    end
    @debug @val lb
    @debug @val ub
    @debug "$([round(e.x,3) for e in collect(pq[t.root])])"
    xn = clip_front(pq[r], line(r, -σ(r)), 0.0)
    return dp_tree_backtrace(xn, t, lb, ub)
end
