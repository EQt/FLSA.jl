### Some algorithms for tree (sub)graphs
import Graph: AbstractGraph

type TreeSubGraph{V,E}
    graph::Graph.AbstractGraph{V,E}
    edges::Vector{E}
    root::V
    parent::Vector{Int}
    children::Vector{Vector{Int}}
    dfs_order::Vector{Int}
    edge_index::Dict{@compat(Tuple{Int,Int}),Int}
end

typealias ITreeSubGraph TreeSubGraph{Int,IEdge}
typealias Tree ITreeSubGraph

"""Tree(n): create an empty tree containing n nodes"""
function Tree(n::Int)
    graph = Graphs.simple_edgelist(n, IEdge[])
    Tree(graph, [], n, [], Array{Array{Int64,1},1}(), [],
         Dict{Tuple{Int64,Int64},Int64}())
end



function Tree(pi::Vector{Int})
    n = length(pi)
    edges = Vector{IEdge}()
    sizehint!(edges, n-1)
    children = [Vector{Int}() for i=1:n]
    root = 0
    k = 1
    for (i,j) in enumerate(pi)
        if i == j
            root = i
        else
            push!(edges, IEdge(k, i,j))
            k += 1
            push!(children[j], i)
        end
    end
    graph = Graph.simple_edgelist(n, edges; is_directed=false)
    dfs_order = Vector{Int}(n)
    stack = Vector{Int}()
    sizehint!(stack, n)
    push!(stack, root)
    dfs_num = 1
    while !isempty(stack)
        v = pop!(stack)
        dfs_order[dfs_num] = v
        dfs_num += 1
        append!(stack, children[v])
    end
    Tree(graph, edges, root, pi, children, dfs_order,
              Dict{Tuple{Int64,Int64},Int64}())
end


"""Create a tree out of a parent array"""
function create_tree(parent::Vector{Int})
    n = length(parent)
    root = 0
    for (i,p) in enumerate(parent)
        if i == p
            root = i
            break
        end
    end
    E = collect(filter(e -> e[1] != e[2], enumerate(parent)))
    children = [Int[] for i=1:n]
    for (i,p) in E
        push!(children[p], i)
    end
    dfs_order = zeros(Int, n)
    stack = [root]
    dfs_num = 1
    while !isempty(stack)
        v = pop!(stack)
        dfs_order[dfs_num] = v
        dfs_num += 1
        append!(stack, children[v])
    end
    E = [(min(s,t), max(s,t)) for (s,t) in E]
    edges = [Edge{Int}(i,min(s,t), max(s,t)) for (i,(s,t)) in enumerate(E)]
    graph = EdgeList{Int,Edge{Int}}(false, collect(1:n), edges)
    edge_index = Dict((e,i) for (i,e) in enumerate(E))
    return Tree(graph, edges, root, parent, children, dfs_order, edge_index)
end
    

"""Compute a subtree with given `root` node by extracting `edges` of `graph`"""
function subtree(graph, edges, root)
    @assert length(edges) > 0
    local n = num_vertices(graph)
    parent = fill(-1, n)
    children = fill(Int[], n)
    dfs_order = fill(-1, n)
    ei = Dict{@compat(Tuple{Int,Int}),Int}()
    tree = TreeSubGraph(graph, edges, root, parent, children, dfs_order, ei)
    init_tree(tree)
    tree
end

"""Initialize all the values by performing a dfs according to the edges"""
function init_tree(t::TreeSubGraph)
    local n = num_vertices(t.graph)
    local neighbors = [Int[] for i=1:n ]
    @inbounds for e in t.edges
        is, it = source(e), target(e)
        push!(neighbors[is], it)
        push!(neighbors[it], is)
    end
    local dfs_num::Int = 1
    iroot = t.root
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
    @inbounds for (k,e) in enumerate(t.edges)
        i, j = source(e), target(e)
        i, j = min(i,j), max(i,j)
        t.edge_index[(i,j)] = k
    end
    t
end


"""Give the indices of an pre-order run through the nodes of `t`, exluding root"""
@inline preorder{V,E}(t::TreeSubGraph{V,E}) = t.dfs_order[2:end]

"""Give the indices of an post-order run through the nodes of `t`, excluding root"""
@inline postorder{V,E}(t::TreeSubGraph{V,E}; include_root::Bool=false) =
    if include_root
        t.dfs_order[end:-1:1]
    else
        t.dfs_order[end:-1:2]
    end


"""Number of nodes in the graph"""
@inline num_nodes{V,E}(t::TreeSubGraph{V,E}) = length(t.children)


"""Return edge tuple (i,j) with i < j"""
@inline E(c,v) = (min(c,v), max(c,v))


"""Index of the edge `e` in `tree`"""
@inline edge_index(tree::Tree, e::@compat(Tuple{Int,Int})) =
    tree.edge_index[E(e[1],e[2])]

function non_tree(D, edges)
    m = size(D)[1]
    d = ones(m)
    for e in edges
        d[e.index] = 0
    end
    E = sparse(collect(1:m), collect(1:m), d, m, m)
    return E * D
end

function tree_part(D, edges)
    m = size(D)[1]
    n = length(edges)
    d = ones(Int, n)
    for (i,e) in enumerate(edges)
        d[i] = e.index
    end
    E = sparse(collect(1:n), d, ones(Int, n), n, m)
    return E * D
end
