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

@inline pix2ind(i::Int, j::Int, g::ImgGraph) = pix2ind(i, j, g.n1)
@inline pix2ind(i::Int, j::Int, n1::Int) = i + (j-1)*n1

@inline ind2pix(i::Int, g::ImgGraph) = nothing


function img_graph(n1::Int, n2::Int, dir = [((1,1), 1.0)])
    print(dir)
    n = n1 * n2
    E = [IEdge(1,n1,n2)]
    m = 0
    for d in dir
        e = d[1]
        m += (n1-e[1])*(n2-e[2]) + (n1-e[2])*(n2-e[1])
    end
    I = zeros(Int, 2m)
    J = zeros(Int, 2m)
    W = zeros(Float64, 2m)
    m = 0
    for d in dir
        e = d[1]
        c = (n1 + n2)*sum(d[1])
        for j = 1:n2-e[2]
            for i = 1:n1-e[1]
                k = m + 2*(i + (j-1)*(n1-e[1])) - 1
                I[k] = k
                J[k] = pix2ind(i,j, n1)
                W[k] = +d[2]
                k += 1
                I[k] = k-1 # same edge
                J[k] = pix2ind(i+e[1], j+e[2], n1)
                W[k] = -d[2]
                @debug @val k
                @debug @val pix2ind(i+e[1], j+e[2], n1)
            end
        end
        m += 2*(n1-e[1])*(n2-e[2])
        for j = 1:n2-e[1]
            for i = 1+e[2]:n1
                k = m + 2*(i - e[2] + (j-1)*(n1-e[2])) - 1
                I[k] = k
                J[k] = pix2ind(i,j, n2)
                W[k] = +d[2]
                k += 1
                I[k] = k-1 # same edge
                J[k] = pix2ind(i-e[2], j+e[1], n2)
                W[k] = -d[2]
            end
        end
        m += 2*(n1-e[2])*(n2-e[1])
    end
    println(I)
    D = sparse(I, J, W)
    G = simple_edgelist(n1*n2, E; is_directed=false)
    ImgGraph(n1, n2, 1.0, G, D)
end
