"""Naive implementation of the dynamic programming algorithm on trees."""


"""Contains all information for the calculation"""
type DPVisitor <: AbstractGraphVisitor
    y::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    df::Vector{PWL}
    function DPVisitor(y::Vector{Float64})
        n = length(y)
        lb = fill(-Inf, n)
        ub = fill(+Inf, n)
        df = [PWL([]) for i in y]
        new(y, lb, ub, df)
    end
end

"""For testing the graph traversal mechanism"""
function Graphs.close_vertex!(vis::DPVisitor, v)
    println("processing $v")
end


"""After computing the bounds [lb, ub] for each node, compute solution x
   by clipping each edge, i.e. backtracing from root to children"""
function backtrace_dp_tree(x_root, iroot, t, n, ub, lb)
    x = zeros(n)
    x[iroot] = x_root
    for v in post_order(t)
        x[v] = clamp(x[t.parent[v]], lb[v], ub[v])
    end
    return x
end

"""Compute x=FLSA(y, lambda) on a a (sub)tree t, naive PWL implementation"""
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


"""Compute the dual solution to x on a tree (sub)graph"""
function dual_tree_rec(y::Vector{Float64}, x::Vector{Float64}, t::ITreeSubGraph)
    local m = length(t.edges)
    alpha = zeros(length(t.edges))
    local iroot = vertex_index(t.root, t.graph)
    dfs_dual_tree(iroot, alpha, t, x,  deepcopy(y))
    return alpha
end


function dual_tree(y::Vector{Float64}, tree)
    α = zeros(length(tree.edges))
    for c in postorder(tree)
        v = tree.parent[c]
        e = t.edge_index[(min(c,v), max(c,v))]
        α[e] = sign(v-c)*(x[c] - y[c])
        y[v] += sign(c-v)*α[e]
    end
    return α
end


function dfs_dual_tree(v, alpha, t, x, y)
    for c in t.children[v]
        dfs_dual_tree(c, alpha, t, x, y)
        ie = t.edge_index[(min(c,v), max(c,v))]
        alpha[ie] = sign(v-c)*(x[c] - y[c])
        y[v] += sign(c-v)*alpha[ie]
    end
end


"""Record what is happening, when a knot of the PWL is hit"""
immutable Event
    x::Float64      # position
    offset::Float64 # delta offset
    slope::Float64  # delta slope
end
