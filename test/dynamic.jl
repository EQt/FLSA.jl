module TestDynamic
using Base.Test
if VERSION >= v"0.7-" using Random end
import FLSA
import FLSA.@debug
import FLSA.Graph: kruskal_minimum_spantree

# debug(msg) = println("DEBUG: $msg")
const digits = 3

@testset "3 line knots" begin
    sol = [2.0; 2.0; 2.5]
    y = [1.0; 2.0; 3.5]
    t = FLSA.create_tree([2,3,3])
    @debug "*"^70
    xl = FLSA.dp_line_naive(y, 1.0)
    @debug "*"^70
    x = FLSA.dp_tree_naive(y, 1.0, t)
    @test x ≈ sol
    @debug "*"^70
    x = FLSA.dp_tree(y, 1.0, t)
    @test x ≈ sol
end

@testset "random(13)" begin
    srand(13)
    lambda = 1.0
    n1, n2 = 43, 31
    g = FLSA.igraph(FLSA.grid_graph(n1, n2))
    v = FLSA.vertices(g)
    n = length(v)
    w = rand(size(g.edges))
    y = round.(10*rand(n1*n2), 1)
    root = 1

    mst, wmst = kruskal_minimum_spantree(g, w)
    t = FLSA.subtree(g, mst, root)
    x = FLSA.dp_tree_naive(y, lambda, t)
    @debug "*"^80
    x2 = FLSA.dp_tree(y, lambda, t)
    @test round.(x, digits) == round.(x2, digits)
end

end # module
