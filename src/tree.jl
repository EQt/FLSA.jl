### Some algorithms for tree (sub)graphs

using Graphs

type TreeSubGraph{V,E}
    graph::AbstractGraph{V,E}
    edges::Vector{E}
    root::V
    parent::Vector{Int}
    children::Vector{Vector{Int}}
    dfs_order::Vector{Int}
    edge_index::Dict{@compat(Tuple{Int,Int}),Int}
end

typealias ITreeSubGraph TreeSubGraph{Int,IEdge}
typealias Tree ITreeSubGraph


"""Compute a subtree with given `root` node by extracting `edges` of `graph`"""
function subtree(graph, edges, root)
    local n = num_vertices(graph)
    parent = fill(-1, n)
    children = fill(Int[], n)
    dfs_order = fill(-1, n)
    ei = Dict{(Int,Int), Int}()
    tree = TreeSubGraph(graph, edges, root, parent, children, dfs_order, ei)
    init_tree(tree)
    tree
end

function init_tree(t::TreeSubGraph)
    local n = num_vertices(t.graph)
    local neighbors = [Int[] for i=1:n ]
    for e in t.edges
        is, it = vertex_index(source(e), t.graph), vertex_index(target(e), t.graph)
        push!(neighbors[is], it)
        push!(neighbors[it], is)
    end
    local dfs_num::Int = 1
    iroot = vertex_index(t.root, t.graph)
    t.parent[iroot] = iroot
    local stack = [iroot]
    while !isempty(stack)
        v = pop!(stack)
        t.dfs_order[dfs_num] = v
        dfs_num += 1
        t.children[v] = neighbors[v]
        for u in t.children[v]
            t.parent[u] = v
            iv = findfirst(neighbors[u], v)
            deleteat!(neighbors[u], iv)
            push!(stack, u)
        end
    end
    for (k,e) in enumerate(t.edges)
        i, j = source(e), target(e)
        i, j = vertex_index(i, t.graph), vertex_index(j, t.graph)
        i, j = min(i,j), max(i,j)
        t.edge_index[(i,j)] = k
    end
    t
end


"""Give the indices of an pre-order run through the nodes of `t`"""
preorder{V,E}(t::TreeSubGraph{V,E}) = t.dfs_order[end:-1:1]


"""Return edge tuple (i,j) with i < j"""
@inline E(c,v) = (min(c,v), max(c,v))


"""Index of the edge `e` in `tree`"""
@inline edge_index(tree::Tree, e::@compat(Tuple{Int,Int})) =
    tree.edge_index[E(e[1],e[2])]
