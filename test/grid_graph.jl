using FLSA
using Base.Test
using Graphs.source, Graphs.target

@test 1 == 1

V23 = [(1,1)
       (1,2)
       (1,3)
       (2,1)
       (2,2)
       (2,3)]

@test FLSA.grid_nodes(2,3) == V23


E23 =[((1,1),(2,1))
      ((1,2),(2,2))
      ((1,3),(2,3))
      ((1,1),(1,2))
      ((1,2),(1,3))
      ((2,1),(2,2))
      ((2,2),(2,3))]

@test FLSA.grid_edges(2,3) == E23

g = FLSA.grid_graph(2,3)

@test g.vertices == V23

for e in g.edges
    @test (source(e), target(e)) in E23
end

@test length(g.edges) == length(E23)
