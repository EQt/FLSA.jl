import Graph: Edge, EdgeList

typealias GridNode @compat Tuple{Int, Int}
typealias GridEdge Graph.Edge{GridNode}
typealias GridGraph Graph.EdgeList{GridNode, GridEdge}

"""Number of rows of the grid"""
num_rows(g::GridGraph) = max(findfirst(x -> x == (1,2), g.vertices)-1, 1)

"""Number of columns of the grid"""
num_columns(g::GridGraph) = div(num_vertices(g), num_rows(g))

"""Vector of nodes of an n1×n2 grid graph"""
function grid_nodes(n1::Integer, n2::Integer)
    V = [(i,j) for i in 1:n1, j in 1:n2]
    reshape(V, length(V))
end

"""Vector of edges of an n1×n2 grid graph"""
function grid_edges(n1::Integer, n2::Integer)
    local Eh = [((i,j), (i,j+1)) for j in 1:n2-1, i in 1:n1]
    local Ev = [((i,j), (i+1,j)) for j in 1:n2,   i in 1:n1-1]
    E = vcat(reshape(Eh, length(Eh)), reshape(Ev, length(Ev)))
end

"""Vector of edges of an n1×n2 diagonal grid graph"""
function dgrid_edges(n1::Integer, n2::Integer)
    local Ep = [((i,j), (i+1,j+1)) for j in 1:n2-1, i in 1:n1-1]
    local En = [((i-1,j+1), (i,j)) for j in 1:n2-1,   i in 2:n1]
    E = vcat(reshape(Ep, length(Ep)), reshape(En, length(En)))
end

"""Construct an incidence list representation of the n1×n2 grid graph"""
function grid_graph(n1::Int, n2::Int)
    V = grid_nodes(n1, n2)
    E = grid_edges(n1, n2)
    E = [GridEdge(i, u, v) for (i, (u, v)) in enumerate(E)]
    edgelist(V, E; is_directed=false)
end

"""Construct an incidence list representation of the n1×n2 diagonal grid graph"""
function dgrid_graph(n1::Int, n2::Int)
    V = grid_nodes(n1, n2)
    E = dgrid_edges(n1, n2)
    E = [GridEdge(i, u, v) for (i, (u, v)) in enumerate(E)]
    edgelist(V, E; is_directed=false)
end
    
