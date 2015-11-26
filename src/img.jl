typealias Pixel Tuple{Int,Int}
import Graphs.IEdge

type ImgGraph
    n1::Int
    n2::Int
    Lip::Float64
    graph::AbstractGraph{Int,IEdge}
    D::AbstractSparseArray{Float64, Int, 2}
end

@inline pix2ind(i::Int, j::Int, g::ImgGraph) = i*(g.n2 - 1) + j

@inline ind2pix(i::Int, g::ImgGraph) = nothing


function img_graph(n1::Int, n2::Int, dir = [(Pixel(1,1), 1.0)])
    n = n1 * n2
    E = [IEdge(1,n1,n2)]
    D = sparse([1.0])
    G = simple_edgelist(n1*n2, E; is_directed=false)
    ImgGraph(n1, n2, 1.0, G, D)
end

