using Graphs

typealias Pixel Tuple{Int,Int}
import Graphs.IEdge

type ImgGraph
    n1::Int
    n2::Int
    Lip::Float64
    graph::AbstractGraph{Int,IEdge}
    D::AbstractSparseArray{Float64, Int, 2}
end

@inline pix2ind(i::Int, j::Int, g::ImgGraph) = pix2ind(i, j, g.n2)
@inline pix2ind(i::Int, j::Int, n2::Int) = i*(n2-1) + j


@inline ind2pix(i::Int, g::ImgGraph) = nothing


function img_graph(n1::Int, n2::Int, dir = [((1,1), 1.0)])
    print(dir)
    n = n1 * n2
    E = [IEdge(1,n1,n2)]
    m = 0
    for d in dir
        e = d[1]
        m += n - n1*e[1]
        m += n - n2*e[2]
    end
    @debug "m = $m"
    I = zeros(Int, 2m)
    J = zeros(Int, 2m)
    W = zeros(Float64, 2m)
    m = 0
    for d in dir
        e = d[1]
        for i = 1:n1-e[1]
            for j = 1:n2-e[2]
                k = m + 2 * pix2ind(i,j, n2)
                I[k] = k
                J[k] = pix2ind(i,j, n2)
                W[k] = +d[2]
                I[k] = k+1
                J[k] = pix2ind(i+e[1], j+e[2], n2)
                W[k] = -d[2]
                # rot 90 degrees
            end
        end
    end
    println(I)
    D = sparse(I, J, W)
    G = simple_edgelist(n1*n2, E; is_directed=false)
    ImgGraph(n1, n2, 1.0, G, D)
end

