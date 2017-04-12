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


"""
    dual_tree(y, tree)

Compute dual solution α on a tree, such that `y = D'*α`.
Due to computational issues, `α` is sorted in parent-order."""
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

function dual_tree0(y::Vector{ℝ}, tree, mu)
    α = zeros(length(tree.edges))
    for c in postorder(tree)
        v = tree.parent[c]
        e = edge_index(tree, (c,v))
        α[e]  = sign(c-v) * mu(c) * y[c]
        y[v] += sign(c-v) * α[e]
    end
    return α
end


dual_tree(y::Vector{Float64}, x::Vector{Float64}, tree) = dual_tree(x-y, tree)
dual_tree0(y::Vector{Float64}, x::Vector{Float64}, tree, mu) =
    dual_tree0(x-y, tree, mu)


"""Represent a line segment"""
type LineSegment
    slope::ℝ
    offset::ℝ
    function LineSegment(s::ℝ, o::ℝ)
        # @assert abs(s) > 1e-15
        new(s,o)
    end
end


@inline event_time(e::Event) = e.x

include("queues.jl")

"""Find x, such that t = slope*x + offset"""
@inline find_x(t::ℝ, slope::ℝ, offset::ℝ, def::ℝ = 0.0) =
    abs(slope) < 1e-15 ? def : (t - offset)/slope

@inline find_x(t::ℝ, ls::LineSegment, def::ℝ = 0.0) = 
    find_x(t, ls.slope, ls.offset, def)


"""Clip PWL represented by `pq` from negative until `t`, starting with `l`"""
function clip_front{Q}(pq::Q, l::LineSegment, t::ℝ)
    x = find_x(t, l, -Inf)
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
    x = find_x(t, l, +Inf)
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

"""FLSA on a tree, computed by a *dynamic programming* algorithm"""
dp_tree(y::Vector{ℝ}, λ::ℝ, t::Tree) =
    dp_tree(y, i->λ, i->1.0, t)
dp_tree(y::Matrix{ℝ}, λ::ℝ, t::Tree) =
    reshape(dp_tree(vec(y), λ, t), size(y))
dp_tree(y::Vector{ℝ}, λ::Vector{ℝ}, t::Tree) =
    dp_tree(y, i->λ[i], i->1.0, t)
dp_tree(y::Vector{ℝ}, λ::Vector{ℝ}, µ::Vector{ℝ}, t::Tree) =
    dp_tree(y, i->λ[i], i->µ[i], t)

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
