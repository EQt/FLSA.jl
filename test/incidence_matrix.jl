n1, n2 = 2, 3
G = grid_graph(n1, n2)
D = incidence_matrix(G)

n = num_vertices(G)
m = num_edges(G)

# test shape
@test ndims(D) == 2
@test size(D, 1) == m
@test size(D, 2) == n

D23 = [# vertices
       # 1  2  3  4  5  6
        +1  0 -1  0  0  0; # edge 1
         0  0 +1  0 -1  0; # edge 2
         0 +1  0 -1  0  0; # edge 3
         0  0  0 +1  0 -1; # edge 4
        +1 -1  0  0  0  0; # edge 5
         0  0 +1 -1  0  0; # edge 6
         0  0  0  0 +1 -1; # edge 7
       ]

@test full(D) == D23
