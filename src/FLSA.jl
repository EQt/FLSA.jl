module FLSA

using Graphs
using Compat: findlast, @compat

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
       inplace

include("debug.jl")
include("flsa.jl")
include("igraph.jl")
include("grid_graph.jl")
include("incidence_matrix.jl")
include("admm.jl")
include("fista.jl")
include("pwl.jl")
include("tree.jl")
include("dynamic.jl")
include("line.jl")
include("deque.jl")
include("depq.jl")
include("inplace.jl")

end # module
