# Legacy import of needed Graphs algorithms

module Graph

abstract AbstractGraph{V, E}

immutable Edge{V}
    index::Int
    source::V
    target::V
end
typealias IEdge Edge{Int}

type GenericEdgeList{V,E,VList,EList} <: AbstractGraph{V,E}
    is_directed::Bool
    vertices::VList
    edges::EList
end

typealias SimpleEdgeList{E} GenericEdgeList{Int,E,UnitRange{Int},Vector{E}}

# construction

simple_edgelist{E}(nv::Integer, edges::Vector{E}; is_directed::Bool=true) =
    SimpleEdgeList{E}(is_directed, intrange(nv), edges)

edgelist{V,E}(vertices::Vector{V}, edges::Vector{E}; is_directed::Bool=true) =
    EdgeList{V,E}(is_directed, vertices, edges)

is_directed(g::GenericEdgeList) = g.is_directed

num_vertices(g::GenericEdgeList) = length(g.vertices)
vertices(g::GenericEdgeList) = g.vertices

num_edges(g::GenericEdgeList) = length(g.edges)
edges(g::GenericEdgeList) = g.edges
edge_index(e, g::GenericEdgeList) = edge_index(e)


end
