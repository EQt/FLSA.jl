using FLSA.grid_graph

n1, n2 = 2, 3
G = grid_graph(n1, n2)
D = incidence_matrix(G)

n = num_vertices(G)
m = num_edges(G)

# test shape
@test ndims(D) == 2
@test size(D, 1) == m
@test size(D, 2) == n
