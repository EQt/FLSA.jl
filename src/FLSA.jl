module FLSA

using Compat: findlast, @compat, @inline
using Graphs
using DataStructures

const ∞ = Inf
typealias ℝ Float64
typealias Node Int
typealias IncMat SparseMatrixCSC{Float64,Int}


export num_vertices,
       igraph,
       grid_graph,
       incidence_matrix,
       vertex_index,
       norm2,
       conjugate_gradient,
       admm,
       fista,
       PWL,
       Knot,
       call,
       find_x,
       clip, clip_x,
       dp_tree,
       dp_tree_naive,
       dual_tree,
       duality_gap,
       min_knot!,
       max_knot!,
       PWLNode,
       PWLTree,
       PDeQue,
       inplace,
       string,
       img_graph,
       mst_tree

"""Record what is happening, when a knot of the PWL is hit"""
immutable Event
    x::Float64      # position
    slope::Float64  # delta slope
    offset::Float64 # delta offset
    function Event(x, s, o)
        @assert isfinite(x)
        @assert abs(s) > 1e-16
        new(x, s, o)
    end
end

include("debug.jl")
include("inplace.jl")
include("deque.jl")
include("flsa.jl")
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
include("mst_tree.jl")

end # module
