"""Some algorithms for tree (sub)graphs"""

using Graphs

type TreeSubGraph{V,E}
    graph::AbstractGraph{V,E}
    edges::Vector{E}
    root::V
    parent::Vector{Int}
    children::Vector{Vector{Int}}
    dfs_order::Vector{Int}
end

function subtree(graph, edges, root)
    local n = num_vertices(graph)
    parent = fill(-1, n)
    children = fill(Int[], n)
    dfs_order = fill(-1, n)
    tree = TreeSubGraph(graph, edges, root, parent, children, dfs_order)
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
    for i=1:n println("neighbors[$i] = $(neighbors[i])") end
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
            println("iv = $iv, length(nei) = $(length(neighbors[u]))")
            deleteat!(neighbors[u], iv)
            push!(stack, u)
        end
    end
    t
end

                           
