# Some tests how to compute MST in julia

module TestTreeFunctions

import FLSA
import Graph
using Base.Test

@testset "construct a tree" begin
    @test FLSA.Tree([1,1,2,2]) != nothing
end

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


@testset "very basic" begin
    c = [0.0, -1.5, 0.5]
    t = FLSA.create_tree([1,1,1])
    a = FLSA.dual_tree(c, t)
    @test a ≈ [1.5, -0.5]
end


@testset "tree_part, dual_tree" begin
    @test norm(img.D * ones(n)) < 1e-12

    D = FLSA.tree_part(img.D, mst)
    @test size(D, 1) == n - 1
    @test size(D, 2) == n
    for i = 1:n-1
        @test nnz(D[i, :]) == 2
    end
    @test norm(D * ones(n)) < 1e-12

    c = rand(size(D, 2))
    c -= mean(c)
    @static if VERSION < v"0.6"
        @test abs(mean(c)) <= 1e-16
    else
        @test mean(c) ≈ 0 atol=1e-16
    end
    a = FLSA.dual_tree(deepcopy(c), tm)
    @test c ≈ D'*a
end


@testset "test ||alpha||∞ <= lambda" begin
    x = FLSA.dp_tree(y, lambda, tm)
    alpha = FLSA.dual_tree(y, x, tm)
    @test maximum(abs.(alpha)) <= lambda + 1e-9
end


end
