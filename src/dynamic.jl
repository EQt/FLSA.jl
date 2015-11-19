import Base.string

"""
After computing the bounds [lb, ub] for each node, compute optimal solution `x`
by clipping each edge, i.e. backtracing from root to children
"""
function backtrace_dp_tree(xr::ℝ, t::Tree, ub, lb)
    x = zeros(num_nodes(t))
    x[t.root] = xr
    for v in postorder(t)
        x[v] = clip(x[t.parent[v]], lb[v], ub[v])
    end
    return x
end


"""Compute x=FLSA(y, λ) on a a (sub)tree t, naive PWL implementation"""
function dp_tree_naive(y::Vector{ℝ}, λ, µ, t::Tree)
    n = length(y)
    df = [ PWL(0.0, -µ(i)*y[i]; slope=µ(i)) for i=1:n ]
    lb, ub = zeros(n), zeros(n)
    for c in preorder(t)
        v = t.parent[c]
        lb[c] = find_x(df[c], -λ[c])
        ub[c] = find_x(df[c], +λ[c])
        df[v] += clip(df[c], lb[c], ub[c])
    end
    xr = find_x(df[t.root], 0.0)
    backtrace_dp_tree(xr, t, lb, ub)
end


"""Compute dual solution α on a tree, such that y = D' α"""
function dual_tree(y::Vector{ℝ}, tree)
    α = zeros(tree.edges)
    for c in postorder(tree)
        v = tree.parent[c]
        e = edge_index(tree, (c,v))
        α[e]  = sign(v-c) * y[c]
        y[v] += sign(c-v) * α[e]
    end
    return α
end


"""Represent a line segment"""
type LineSegment
    slope::ℝ
    offset::ℝ
end


"""Record what is happening, when a knot of the PWL is hit"""
immutable Event
    x::Float64      # position
    slope::Float64  # delta slope
    offset::Float64 # delta offset
    # Event(x, l::LineSegment) = new(x, l.slope, l.offset)
end

"""Print them more readable"""
string(e::Event) = "$(e.x) @ $(e.slope)x + $(e.slope))"


"""Find x, such that t = slope*x + offset"""
@inline find_x(t::ℝ, slope::ℝ, offset::ℝ) = (t - offset)/slope

@inline find_x(t::ℝ, ls::LineSegment) = find_x(t, ls.slope, ls.offset)
