# Legacy import of needed Graphs algorithms

module Graph

export edgelist,
    is_directed,
    num_vertices,
    num_edges,
    source,
    target,
    implements_edge_list,
    implements_adjacency_list,
    edges,
    vertices,
    kruskal_minimum_spantree,
    simple_edgelist

import Base.isless
using DataStructures: IntDisjointSets, in_same_set, union!, num_groups

abstract type AbstractGraph{V,E} end

immutable Edge{V}
    index::Int
    source::V
    target::V
end
const IEdge = Edge{Int}

import Base.show
function show(io::IO, e::Edge)
    print(io, "edge [$(e.index)]: $(e.source) -- $(e.target)")
end

source(e::Edge) = e.source
target(e::Edge) = e.target
source{V}(e::Edge{V}, g::AbstractGraph{V}) = e.source
target{V}(e::Edge{V}, g::AbstractGraph{V}) = e.target

type GenericEdgeList{V,E,VList,EList} <: AbstractGraph{V,E}
    is_directed::Bool
    vertices::VList
    edges::EList
end

SimpleEdgeList{E} = GenericEdgeList{Int,E,UnitRange{Int},Vector{E}}
EdgeList{V,E} = GenericEdgeList{V,E,Vector{V},Vector{E}}

# construction
intrange(n::Integer) = 1:convert(Int,n)

simple_edgelist{E}(nv::Integer, edges::Vector{E}; is_directed::Bool=true) =
    SimpleEdgeList{E}(is_directed, intrange(nv), edges)

edgelist{V,E}(vertices::Vector{V}, edges::Vector{E}; is_directed::Bool=true) =
    EdgeList{V,E}(is_directed, vertices, edges)

is_directed(g::GenericEdgeList) = g.is_directed

num_vertices(g::GenericEdgeList) = length(g.vertices)
vertices(g::GenericEdgeList) = g.vertices

num_edges(g::GenericEdgeList) = length(g.edges)
edges(g::GenericEdgeList) = g.edges
@inline edge_index(e::Edge) = e.index
@inline edge_index(e, g::GenericEdgeList) = edge_index(e)

@inline vertex_index(v::Int) = v
@inline vertex_index(v, g::GenericEdgeList) = vertex_index(v)

implements_adjacency_list(g::GenericEdgeList) = false
implements_edge_list(g::GenericEdgeList) = true

# Kruskal's algorithm for minimum spanning tree/forest

immutable WeightedEdge{E,W}
    edge::E
    weight::W
end
isless{E,W}(a::WeightedEdge{E,W}, b::WeightedEdge{E,W}) = a.weight < b.weight

abstract type AbstractEdgePropertyInspector{T} end

type VectorEdgePropertyInspector{T} <: AbstractEdgePropertyInspector{T}
  values::Vector{T}
end


edge_property{T,V}(visitor::VectorEdgePropertyInspector{T},
                   e, g::AbstractGraph{V}) =
    visitor.values[edge_index(e, g)]


function collect_weighted_edges{V,E,W}(graph::AbstractGraph{V,E}, weights::AbstractEdgePropertyInspector{W})

    # edge_property_requirement(weights, graph)

    wedges = Array(WeightedEdge{E,W}, 0)
    sizehint!(wedges, num_edges(graph))

    if implements_edge_list(graph)
        for e in edges(graph)
            w = edge_property(weights, e, graph)
            push!(wedges, WeightedEdge(e, w))
        end

    elseif implements_vertex_list(graph) && implements_incidence_list(graph)
        for v in vertices(graph)
            for e in out_edges(v, graph)
                w = edge_property(weights, e, graph)
                push!(wedges, WeightedEdge(e, w))
            end
        end
    else
        throw(ArgumentError("graph must implement either edge_list or incidence_list."))
    end

    return wedges
end


# select edges from a sorted list of weighted edges
function kruskal_select{V,E,W}(
    graph::AbstractGraph{V,E},
    sorted_wedges::AbstractVector{WeightedEdge{E,W}},
    K::Integer)

    n = num_vertices(graph)
    re = Array(E, 0)
    rw = Array(W, 0)

    if n > 1
        dsets = IntDisjointSets(n)
        sizehint!(re, n-1)
        sizehint!(rw, n-1)

        ui::Int = 0
        vi::Int = 0

        for we in sorted_wedges
            e::E = we.edge
            ui = vertex_index(source(e, graph), graph)
            vi = vertex_index(target(e, graph), graph)

            if !in_same_set(dsets, ui, vi)
                union!(dsets, ui, vi)
                push!(re, e)
                push!(rw, we.weight)
            end

            if num_groups(dsets) <= K
                break
            end
        end
    end

    return (re, rw)
end

function kruskal_minimum_spantree(graph::AbstractGraph, eweights::AbstractEdgePropertyInspector; K::Integer=1)

    # collect & sort edges

    wedges = collect_weighted_edges(graph, eweights)
    sort!(wedges)

    # select the tree edges
    kruskal_select(graph, wedges, K)
end


function kruskal_minimum_spantree(graph::AbstractGraph, eweights::AbstractVector; K::Integer=1)
    visitor::AbstractEdgePropertyInspector = VectorEdgePropertyInspector(eweights)
    kruskal_minimum_spantree(graph, visitor, K=K)
end

end
