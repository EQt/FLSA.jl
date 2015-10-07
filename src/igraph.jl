"""Integer Graphs to avoid slow vertex_index calls"""

import Graphs.IEdge

typealias IGraph GenericEdgeList{Int,IEdge,Vector{Int},Vector{IEdge}}

function igraph{V,E}(g::AbstractGraph{V,E},
                     vmap::Dict{V,Int} = Dict([(v, i) for (i,v) in enumerate(vertices(g))]))
    ed = [IEdge(e.index, vmap[source(e)], vmap[target(e)]) for e in edges(g)]
    for (i,e) in enumerate(ed)
        u,v = source(e), target(e)
        if (u > v)
            ed[i] = IEdge(e.index, v, u)
        end
    end
    edgelist(collect(1:num_vertices(g)), ed; is_directed=false)
end

if !isdefined(current_module(), :vertex_index)
    # avoid warning; see answer
    # https://github.com/JuliaLang/julia/issues/7860#issuecomment-51340412

    """Fast implementation of vertex_index"""
    Graphs.vertex_index(v::Int, g::IGraph) = v
end
