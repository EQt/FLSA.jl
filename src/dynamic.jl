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

function Graphs.close_vertex!(vis::DPVisitor, v)
    println("processing $v")
end



"""After computing the bounds [lb, ub] for each node, compute solution x
   by clipping each edge, i.e. backtracing from root to children"""
function backtrace_dp_tree(vis, iroot, t, y)
    x_root = find_x(vis.df[iroot], 0)
    x = zeros(y)
    x[iroot] = x_root
    for v in t.dfs_order[2:end]
        x[v] = clamp(x[t.parent[v]], vis.lb[v], vis.ub[v])
    end
    return x
end

"""Compute x=FLSA(y, lambda) on a a (sub)tree t, naive PWL implementation"""
function dp_tree_naive(y::Vector{Float64}, lambda::Float64, t::TreeSubGraph,
                       vis::DPVisitor = DPVisitor(y))
    for i in 1:num_vertices(t.graph)
        vis.df[i] = PWL([Knot(-INF, -INF - vis.y[i]),
                         Knot(+INF, +INF - vis.y[i])])
    end
    for c in t.dfs_order[end:-1:1]
        v = t.parent[c]
        vis.lb[c] = find_x(vis.df[c], -lambda)
        vis.ub[c] = find_x(vis.df[c], +lambda)
        vis.df[c] = clip(vis.df[c], vis.lb[c], vis.ub[c])
        vis.df[v] += vis.df[c]
    end

    local iroot = vertex_index(t.root, t.graph)
    backtrace_dp_tree(vis, iroot, t, y)
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
function dual_tree(y::Vector{Float64}, x::Vector{Float64}, t::ITreeSubGraph)
    local m = length(t.edges)
    alpha = zeros(length(t.edges))
    local iroot = vertex_index(t.root, t.graph)
    dfs_dual_tree(iroot, alpha, t, x,  deepcopy(y))
    return alpha
end

function dfs_dual_tree(v, alpha, t, x, y)
    for c in t.children[v]
        dfs_dual_tree(c, alpha, t, x, y)
        ie = t.edge_index[(min(c,v), max(c,v))]
        alpha[ie] = sign(v-c)*(x[c] - y[c])
        y[v] += sign(c-v)*alpha[ie]
    end
end
