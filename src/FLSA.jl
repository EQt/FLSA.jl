# __precompile__()

if VERSION >= v"0.7-"
    using Base.include
else
    include("Graphs.jl")
end


module FLSA

# Fix https://github.com/JuliaLang/DataStructures.jl/issues/200
# if Pkg.installed("DataStructures") == v"0.4.3" && VERSION.minor == 5
#     Pkg.checkout("DataStructures", "master")
# end

const CUSTOM_PRINTER = false

using Compat: findlast, @compat, @inline, String
if VERSION >= v"0.7-"
    include("Graphs.jl")
    using .Graph
    using .Graph: IEdge
    using Printf
else
    using Graph
    using Graph: IEdge
end

using DataStructures

const ∞ = Inf
const ℝ = Float64
const Node = Int
const IncMat = SparseMatrixCSC{Float64,Int}

export num_vertices,
       grid_graph,
       incidence_matrix,
       norm2,
       conjugate_gradient,
       admm,
       fista,
       dp_tree,
       dp_tree_naive,
       dual_tree,
       duality_gap,
       string,
       img_graph,
       mst_tree

@compat abstract type Element end

"""Record what is happening, when a knot of the PWL is hit"""
struct Event <: Element
    x::Float64      # position
    slope::Float64  # delta slope
    offset::Float64 # delta offset
    function Event(x, s, o)
        # @assert isfinite(x)
        # @assert abs(s) > 1e-16
        new(x, s, o)
    end
end

import Base.isless
isless(e1::Event, e2::Event) = isless(e1.x, e2.x)


if CUSTOM_PRINTER
    import Base.string
    function string(e::Event)
        @sprintf "@%f : Δs = %f, Δo = %f" e.x e.slope e.offset
    end


    import Base.show
    function show(io::IO, e::Event)
        print(io, string(e))
    end
end

include("debug.jl")
include("deque.jl")
include("utils.jl")
include("igraph.jl")
include("grid_graph.jl")
include("incidence_matrix.jl")
include("tree.jl")
include("img.jl")
include("admm.jl")
include("fista.jl")
include("pwl.jl")
include("dynamic.jl")
include("line.jl")
include("mgt.jl")
include("hdf5.jl")

mst_tree = kruskal_minimum_spantree

precompile(img_graph, (Int64,Int64,Array{Tuple{Tuple{Int64,Int64},Float64},1}))
precompile(img_graph, (Int64, Int64, Int64, Float64))
precompile(max_gap_tree, (Array{Float64,1}, ImgGraph))
precompile(fista, (Array{Float64,1},SparseMatrixCSC{Float64,Int64}))
precompile(admm, (Array{Float64,1},SparseMatrixCSC{Float64,Int64}))

end # module
