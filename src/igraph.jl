"""Integer Graphs to avoid slow vertex_index calls"""

typealias IEdge (Int,Int)
typealias IGraph GenericEdgeList{Int,IEdge,Vector{Int},Vector{IEdge}}

function igraph{V,E}(g::AbstractGraph{V,E},
                     vmap::Dict{V,Int} = Dict([(v, i) for (i,v) in enumerate(vertices(g))]))
    ed = [(vmap[source(e)], vmap[target(e)]) for e in edges(g)]
    for (i,(u,v)) in enumerate(ed)
        if (u > v)
            ed[i] = (v, u)
        end
    end
    simple_edgelist(num_vertices(g), ed; is_directed=false)
end

"""Fast implementation of vertex_index"""
Graphs.vertex_index(v::Int, g::IGraph) = v
