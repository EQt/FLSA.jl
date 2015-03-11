"""Vector of nodes of an n1×n2 grid graph"""
function grid_nodes(n1::Integer, n2::Integer)
    V = [(i,j) for j in 1:n2, i in 1:n1]
    reshape(V, length(V))
end

"""Vector of edges of an n1×n2 grid graph"""
function grid_edges(n1::Integer, n2::Integer)
    Eh = [((i,j), (i+1,j)) for j in 1:n2,   i in 1:n1-1]
    Ev = [((i,j), (i,j+1)) for j in 1:n2-1, i in 1:n1]
    E = vcat(reshape(Eh, length(Eh)), reshape(Ev, length(Ev)))
end

typealias GridNode (Int, Int)
typealias GridEdge Graphs.Edge{GridNode}
typealias GridGraph IncidenceList{GridNode, GridEdge}

"""Construct an incidence list representation of the n1×n2 grid graph"""
function grid_graph(n1::Int, n2::Int)
    V = grid_nodes(n1, n2)
    E = grid_edges(n1, n2)
    E = [GridEdge(i, u, v) for (i, (u, v)) in enumerate(E)]
end
