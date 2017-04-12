# Some tests how to compute MST in julia

module TestTreeFunctions

import FLSA
import Graph
using Base.Test

srand(42)

n1, n2 = 3, 2
img = FLSA.img_graph(n1, n2)
g = img.graph
v = Graph.vertices(g)
n = length(v)
y = rand(n1*n2)
mst, wmst = Graph.kruskal_minimum_spantree(g, rand(size(g.edges)))
lambda = 0.1
tm = FLSA.subtree(g, mst, 1)


function test1()
    x = FLSA.dp_tree(y, lambda, tm)
    alpha = FLSA.dual_tree(y, x, tm)
    @test all(abs(alpha) .<= lambda + 1e-9)
end


function test2()
    @test norm(img.D * ones(n)) < 1e-12

    D = FLSA.tree_part(img.D, mst)
    @test size(D, 1) == n - 1
    @test size(D, 2) == n
    for i = 1:n-1
        @test nnz(D[i, :]) == 2
    end
    @test norm(D * ones(n)) < 1e-12
end

end

TestTreeFunctions.test1()
TestTreeFunctions.test2()
