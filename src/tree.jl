"""Some algorithms for tree (sub)graphs"""

using Graphs

type TreeSubGraph{V,E}
    g::AbstractGraph{V,E}
    edges::Vector{E}
    root::V
    parent::Vector{Int}
    children::Vector{Vsdfasdfaector{Int}}
    dfs_order::Vector{Int}
end

function subtree(graph, edges, root)
    local n = num_vertices(g)
    local neighbors = fill(Int[], n)
    parent = fill(-1, n)
    children = fill(Int[], n)
    dfs_order = fill(-1, n)
    for e in edges
        is, it = vertex_index(source(e), graph), vertex_index(target(e), graph)
        push!(neighbors[is], it)
        push!(neighbors[it], is)
    end
    local dfs_num::Int = 1
    local stack = [vertex_index(root, graph)]
    while !isempty(stack)
        v = pop!(stack)
        dfs_order[dfs_num] = v
        dfs_num += 1
        children[v] = neighbors[v]
        for u in children[v]
            parent[u] = v
            push!(stack, u)
        end
    end
    TreeSubGraph(graph, edges, root, parent, children)
end
