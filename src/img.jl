using Graphs
import Graphs.IEdge
import Base.norm

typealias Pixel Tuple{Int,Int}

norm2(x::Pixel) = x[1]^2 + x[2]^2

function __init__()
    norm(x::Pixel) = sqrt(norm2(x))
end


type ImgGraph
    n1::Int
    n2::Int
    Lip::Float64
    lambda::Vector{Float64}
    graph::IGraph
    D::IncMat
end

# overload for convinience
gap_vec(y, alpha, grid::ImgGraph) = gap_vec(y, alpha, grid.D)
duality_gap(y, alpha, grid::ImgGraph) = sum(gap_vec(y, alpha, grid))
dp_tree(y::Vector{ℝ}, g::ImgGraph, t::Tree) = dp_tree(y, g.lambda, t)
dp_tree(y::Matrix{ℝ}, g::ImgGraph, t::Tree) =
    reshape(dp_tree(y[:], g.lambda, t), g.n1, g.n2)
fista(y::Matrix{ℝ}, g::ImgGraph; ps...) =
    reshape(fista(y[:], g.D; L=g.Lip, ps...), g.n1, g.n2)

@inline pix2ind(i::Int, j::Int, g::ImgGraph) = pix2ind(i, j, g.n1)
@inline pix2ind(i::Int, j::Int, n1::Int) = i + (j-1)*n1

@inline ind2pix(i::Int, g::ImgGraph) = nothing


"""
Create an grid graph for an images of `n1`*`n2` pixels.

`dn`  is the number of directions to go on from an internal pixel node

`lam` is a scaling factor.
"""
function img_graph(n1::Int, n2::Int, dn::Int, lam::Float64)
    println("n1 = $n1, n2 = $n2")
    g = img_graph(n1, n2, dn)
    g.D *= lam
    g.Lip *= lam^2
    g
end


"""
`ds` species the directions explicitly (vector d is scaled by 1/norm(d))
"""
img_graph(n1::Int, n2::Int, ds::Vector{Pixel}) =
    img_graph(n1, n2, [(d, 1/norm(d)) for d in ds])


function img_graph(n1::Int, n2::Int, dn::Int)
    if dn == 1
        img_graph(n1, n2, [(1,0)])
    elseif dn == 2
        img_graph(n1, n2, [(1,0), (1,1)])
    elseif dn == 3
        img_graph(n1, n2, [(1,0), (1,1), (2,0)])
    elseif dn == 4
        img_graph(n1, n2, [(1,0), (1,1), (2, 0), (2,1), (1,2)])
    else
        throw(ArgumentError("dn >= 4 not supported"))
    end
end


function img_graph(n1::Int, n2::Int,
                   dir::Vector{Tuple{Pixel,Float64}}=[((1,0), 1.0)])
    n = n1 * n2
    m = 0
    for d in dir
        e = d[1]
        m += (n1-e[1])*(n2-e[2]) + (n1-e[2])*(n2-e[1])
    end
    @assert m > 0
    I = zeros(Int, 2m)
    J = zeros(Int, 2m)
    W = zeros(Float64, 2m)
    E = [IEdge(0,0,0) for e=1:m]
    lam = zeros(Float64, m)
    m = 0
    for d in dir
        e = d[1]
        @inbounds for j = 1:n2-e[2]
            for i = 1:n1-e[1]
                l = m + i + (j-1)*(n1-e[1])
                k = 2l - 1
                v1 = pix2ind(i,j, n1)
                v2 = pix2ind(i+e[1], j+e[2], n1)
                E[l] = IEdge(l, v1, v2)
                lam[l] = d[2]
                I[k] = l
                J[k] = v1
                W[k] = +d[2]
                k += 1
                I[k] = l # same edge
                J[k] = v2
                W[k] = -d[2]
            end
        end
        m += (n1-e[1])*(n2-e[2])
        @inbounds for j = 1:n2-e[1]
            for i = 1+e[2]:n1
                l = m + i - e[2] + (j-1)*(n1-e[2])
                k = 2l -1
                v1 = pix2ind(i,j, n1)
                v2 = pix2ind(i-e[2], j+e[1], n1)
                E[l] = IEdge(l, v1, v2)
                lam[l] = d[2]
                I[k] = l
                J[k] = v1
                W[k] = +d[2]
                k += 1
                I[k] = l # same edge
                J[k] = v2
                W[k] = -d[2]
            end
        end
        m += (n1-e[2])*(n2-e[1])
    end
    @debug @val I
    @debug @val J
    D = sparse(I, J, W, m, n)
    G = simple_edgelist(n1*n2, E; is_directed=false)
    lmax = maximum(lam)
    Lip = lmax * 2 * 4 * sum([l for (d,l) in dir])
    ImgGraph(n1, n2, Lip, lam, G, D)
end
