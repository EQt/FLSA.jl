"""Vector of nodes of an n1×n2 grid graph"""
function grid_nodes(n1::Integer, n2::Integer)
    V = [(i,j) for j in 1:n2, i in 1:n1]
    reshape(V, length(V))
end

"""Vector of edges of an n1×n2 grid graph"""
function grid_edges(n1::Integer, n2::Integer)
    Eh = [((i,j), (i+1,j)) for i in 1:n1-1, j in 1:n2]
    Ev = [((i,j), (i,j+1)) for i in 1:n1, j in 1:n2-1]
    E = vcat(reshape(Eh, length(Eh)), reshape(Ev, length(Ev)))
end
