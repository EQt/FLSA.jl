module TestMST
include(joinpath("..", "src", "mst.jl"))

import FLSA
using Base.Test

@testset "TestMST" begin
@testset "UnionFind" begin
    n = 5
    uf = UnionFind(n)
    f3, f2 = find(uf, 3), find(uf, 2)
    @test f3 != f2
    unite!(uf, f3, f2)
    f3, f2 = find(uf, 3), find(uf, 2)
    @test f3 == f2
    f1 = find(uf, 1)
    @test f3 != f1
    unite!(uf, f3, f1)
    f3, f1 = find(uf, 3), find(uf, 1)
    @test f3 == f1
    f4, f5 = find(uf, 4), find(uf, 5)
    @test f4 != f5
    @test f4 != f1
    @test f5 != f1
    unite!(uf, f1, f4)
    unite!(uf, f1, f5)
    fi = [find(uf, i) for i in 1:n]
    @test unique(f1) == [f1[1]]
end


@static if !method_exists(sort, (Int, Int))
    Base.sort(a::Int,b::Int) = min(a,b), max(a,b)
end


@testset "ambiguous mst" begin
    n1, n2 = 5, 2
    dn = 2
    lam = 0.1
    root_node = 1
    g = FLSA.img_graph(n1, n2, dn, lam)
    γ = [1.05406e-17,
         0.000666963,
         3.46945e-18,
         0.0,
         0.0,
         3.46945e-18,
         2.1684e-19,
         0.0,
         1.23736e-17,
         0.0,
         2.08167e-17,
         4.33681e-19,
         4.33681e-19,
         1.38778e-17,
         7.11237e-17,
         4.16334e-17,
         4.33681e-19,
         0.0,
         2.1684e-19,
         3.25261e-19,
         4.33681e-19]
    mst, wmst = FLSA.kruskal_minimum_spantree(g.graph, -γ)
    t = FLSA.subtree(g.graph, mst, root_node)
    mst_weight = sum(wmst)
    parent = minimum_spantree(g, -γ, root_node)
    edge_weights = Dict((e.source, e.target) => w
                        for (e,w) in zip(g.graph.edges, γ))
    edge_weights[(root_node, root_node)] = 0.0
    mst_weight2 = -sum(edge_weights[sort(i, p)] for (i,p) in enumerate(parent))
    @test mst_weight ≈ mst_weight2
    @test t.parent != parent

    mst_selected = kruskal_mst(n1*n2, g.graph.edges, γ)
    mst3 = g.graph.edges[mst_selected]
    @test mst == mst3
end


#=
@testset "compare Kruskal and Prim" begin
    # y = TreeInstance.generate(3000)[1]
    # n1, n2 = size(y)
    n1, n2 = 300, 120
    dn = 2
    g = FLSA.img_graph(n1, n2, dn, 0.1)
    m = length(g.graph.edges)
    γ = randn(srand(42), m)
    edges = g.graph.edges

    @time mst_selected = kruskal_mst(n1*n2, edges, γ)
    # @time mst_selected = kruskal_mst(n1*n2, edges, γ)
    @time neigh = compute_undirected_index(mst_selected, edges, n1*n2)
    # @time neigh = compute_undirected_index(view(edges, mst_selected), n1*n2)
    γ .= -γ
    # @time mst2 = minimum_spantree(g, γ)
    @time mst2 = minimum_spantree(g, γ)
end
=#

end
end
