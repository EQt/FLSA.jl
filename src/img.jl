import Graph.IEdge
import Base: norm, convert

struct Pixel
    x::Int
    y::Int
end

convert(::Type{Pixel}, p::Tuple{Int,Int}) = Pixel(p...)
convert(::Type{Tuple{Int,Int}}, p::Pixel) = (p.x, p.y)

norm2(p::Pixel) = p.x^2 + p.y^2
norm(x::Pixel) = sqrt(norm2(x))

# function getindex(p::Pixel, i::Int)
#     if i == 1
#         return p.x
#     elseif i == 2
#         return p.y
#     else
#         throw(ArgumentErrow("i = $i"))
#     end
# end


mutable struct ImgGraph
    n1::Int
    n2::Int
    Lip::Float64
    lambda::Vector{Float64}
    graph::IGraph
    D::IncMat
    dir::Vector{Tuple{Pixel,Float64}}
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
function img_graph(n1::Int, n2::Int, dn::Int, lam::Float64; verbose=true)
    if verbose
        println("n1 = $n1, n2 = $n2")
    end
    g = img_graph(n1, n2, dn)
    g.D *= lam
    g.lambda *= lam
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
        img_graph(n1, n2, Pixel[(1,0)])
    elseif dn == 2
        img_graph(n1, n2, Pixel[(1,0), (1,1)])
    elseif dn == 3
        img_graph(n1, n2, Pixel[(1,0), (1,1), (2,0)])
    elseif dn == 4
        img_graph(n1, n2, Pixel[(1,0), (1,1), (2, 0), (2,1), (1,2)])
    else
        throw(ArgumentError("dn >= 4 not supported"))
    end
end


function img_graph(n1::Int, n2::Int,
                   dir::Vector{Tuple{Pixel,Float64}}=[(Pixel(1,0), 1.0)])
    n = n1 * n2
    m = 0
    for d in dir
        e = d[1]
        m += (n1-e.x)*(n2-e.y) + (n1-e.y)*(n2-e.x)
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
        @inbounds for j = 1:n2-e.y
            for i = 1:n1-e.x
                l = m + i + (j-1)*(n1-e.x)
                k = 2l - 1
                v1 = pix2ind(i,j, n1)
                v2 = pix2ind(i+e.x, j+e.y, n1)
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
        m += (n1-e.x)*(n2-e.y)
        @inbounds for j = 1:n2-e.x
            for i = 1+e.y:n1
                l = m + i - e.y + (j-1)*(n1-e.y)
                k = 2l -1
                v1 = pix2ind(i,j, n1)
                v2 = pix2ind(i-e.y, j+e.x, n1)
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
        m += (n1-e.y)*(n2-e.x)
    end
    @debug @val I
    @debug @val J
    D = sparse(I, J, W, m, n)
    G = simple_edgelist(n1*n2, E; is_directed=false)
    lmax = maximum(lam)
    Lip = lmax * 2 * 4 * sum([l for (d,l) in dir])
    ImgGraph(n1, n2, Lip, lam, G, D, dir)
end
