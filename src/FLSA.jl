module FLSA
  using Graphs

  export num_vertices, grid_graph, incidence_matrix

  include("grid_graph.jl")
  include("incidence_matrix.jl")
end # module
