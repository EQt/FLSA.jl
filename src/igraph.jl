### Integer Graphs to avoid slow vertex_indexo calls

import Graph: IEdge, SimpleEdgeList

const IGraph = Graph.SimpleEdgeList{IEdge}

function igraph{V,E}(g::Graph.AbstractGraph{V,E},
                     vmap::Dict{V,Int} =
                     Dict([(v, i) for (i,v) in enumerate(vertices(g))]))
    ed = [IEdge(e.index, vmap[source(e)], vmap[target(e)]) for e in edges(g)]
    for (i,e) in enumerate(ed)
        u,v = source(e), target(e)
        if (u > v)
            ed[i] = IEdge(e.index, v, u)
        end
    end
    edgelist(collect(1:num_vertices(g)), ed; is_directed=false)
end

