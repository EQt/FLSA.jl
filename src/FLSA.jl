module FLSA
  using Graphs

  export num_vertices,
         grid_graph,
         incidence_matrix,
         norm2,
         conjugate_gradient,
         admm,
         fista,
         PWL,
         Knot,
         call

  include("flsa.jl")
  include("grid_graph.jl")
  include("incidence_matrix.jl")
  include("admm.jl")
  include("fista.jl")
  include("pwl.jl")
end # module
