"""Some algorithms for tree (sub)graphs"""

using Graphs

type TreeSubGraph{V,E}
    g::AbstractGraph{V,E}
    edges::Vector{E}
    root::V
    parent::Vector{Int64}
    children::Vector{Vector{Int64}}
end

function subtree(graph, edges, root)
    local n = num_vertices(g)
    local neighbors = fill(Int64[], n)
    parent = fill(-1, n)
    children = fill(Int64[], n)
    
    TreeSubGraph(graph, edges, root, parent, children)
end
