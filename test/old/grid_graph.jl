@testset "grid_graph" begin

V23 = reshape([(1,1) (1,2) (1,3);
               (2,1) (2,2) (2,3)], 6)

@test FLSA.grid_nodes(2,3) == V23

E23 =[((1,1),(1,2))
      ((1,2),(1,3))
      ((2,1),(2,2))
      ((2,2),(2,3))
      ((1,1),(2,1))
      ((1,2),(2,2))
      ((1,3),(2,3))]

@test FLSA.grid_edges(2,3) == E23

g = FLSA.grid_graph(2,3)

@test is_directed(g) == false
@test g.vertices == V23
@test FLSA.num_columns(g) == 3
@test FLSA.num_rows(g) == 2

@test FLSA.vertex_index(V23[1], g) == 1
@test FLSA.vertex_index(V23[2], g) == 2
@test FLSA.vertex_index(V23[3], g) == 3
@test FLSA.vertex_index(V23[4], g) == 4
@test FLSA.vertex_index(V23[6], g) == 6

for e in g.edges
    @test (source(e), target(e)) in E23
end

@test length(g.edges) == length(E23)

F23 = [((1,1),(2,2))
       ((1,2),(2,3))
       ((1,2),(2,1))
       ((1,3),(2,2))]
@test FLSA.dgrid_edges(2,3) == F23

g = FLSA.dgrid_graph(2,3)
@test g.vertices == V23

end
