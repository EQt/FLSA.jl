## Old code of the recursive implementations

"""Compute the dual solution to x on a tree (sub)graph"""
function dual_tree_rec(y::Vector{Float64}, x::Vector{Float64}, t::ITreeSubGraph)
    local m = length(t.edges)
    alpha = zeros(length(t.edges))
    local iroot = vertex_index(t.root, t.graph)
    dfs_dual_tree(iroot, alpha, t, x,  deepcopy(y))
    return alpha
end


"""Compute x=FLSA(y, lambda) on a a (sub)tree t, recursively"""
function dp_tree_rec(y::Vector{Float64}, lambda::Float64, t::TreeSubGraph)
    vis = DPVisitor(y)
    local iroot = vertex_index(t.root, t.graph)
    dfs_dp_tree(iroot, vis, lambda, t)
    backtrace_dp_tree(vis, iroot, t, y)
end
