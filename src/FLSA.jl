__precompile__()

module FLSA

# Fix https://github.com/JuliaLang/DataStructures.jl/issues/200
if Pkg.installed("DataStructures") == v"0.4.3" && VERSION.minor == 5
    Pkg.checkout("DataStructures", "master")
end

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

end # module
