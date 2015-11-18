"""
After computing the bounds [lb, ub] for each node, compute optimal solution `x`
by clipping each edge, i.e. backtracing from root to children
"""
function backtrace_dp_tree(x_root, iroot, t, n, ub, lb)
    x = zeros(n)
    x[iroot] = x_root
    for v in post_order(t)
        x[v] = clip(x[t.parent[v]], lb[v], ub[v])
    end
    return x
end

"""Compute x=FLSA(y, λ) on a a (sub)tree t, naive PWL implementation"""
function dp_tree_naive{V,E}(y::Vector{Float64}, λ::Float64, t::TreeSubGraph{V,E})
    n = length(y)
    df = [PWL(0, -y[i]; slope=1.0) for i=1:n]
    lb, ub = zeros(n), zeros(n)
    for c in pre_order(t)
        v = t.parent[c]
        lb[c] = find_x(df[c], -λ)
        ub[c] = find_x(df[c], +λ)
        df[c] = clip(df[c], lb[c], ub[c])
        df[v] += df[c]
    end
    local iroot = vertex_index(t.root, t.graph)
    x_root = find_x(vis.df[iroot], 0)
    backtrace_dp_tree(x_root, iroot, t, n, lb, ub)
end


"""Compute x=FLSA(y, lambda) on a a (sub)tree t, recursively"""
function dp_tree_rec(y::Vector{Float64}, lambda::Float64, t::TreeSubGraph)
    vis = DPVisitor(y)
    local iroot = vertex_index(t.root, t.graph)
    dfs_dp_tree(iroot, vis, lambda, t)
    backtrace_dp_tree(vis, iroot, t, y)
end

function dfs_dp_tree(v, vis, lambda, t)
    vis.df[v] = PWL([Knot(-INF, -INF - vis.y[v]),
                     Knot(+INF, +INF - vis.y[v])])
    for c in t.children[v]
        dfs_dp_tree(c, vis, lambda, t)
        vis.lb[c] = find_x(vis.df[c], -lambda)
        vis.ub[c] = find_x(vis.df[c], +lambda)
        vis.df[c] = clip(vis.df[c], vis.lb[c], vis.ub[c])
        vis.df[v] += vis.df[c]
    end
end


"""Return edge tuple (i,j) with i < j"""
@inline E(c,v) = (min(c,v), max(c,v))


"""Compute dual solution α on a tree, such that y = D' α"""
function dual_tree(y::Vector{ℝ}, tree)
    α = zeros(length(tree.edges))
    for c in postorder(tree)
        v = tree.parent[c]
        e = tree.edge_index[E(c,v)]
        α[e]  = sign(v-c) * y[c]
        y[v] += sign(c-v) * α[e]
    end
    return α
end


"""Record what is happening, when a knot of the PWL is hit"""
immutable Event
    x::Float64      # position
    offset::Float64 # delta offset
    slope::Float64  # delta slope
end
