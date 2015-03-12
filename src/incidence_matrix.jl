"""Compute the oriented incidence matrix of g and store in a sparse matrix of type T"""
function incidence_matrix{T<:Number}(g::AbstractGraph, ::Type{T} = Int)
    @graph_requires g vertex_list
    is_directed(g) || error("g must be undirected.")

    local n = num_vertices(g)
    local m = num_edges(g)
    local nnz = 2*m
    I = Array(Int, nnz)
    J = Array(Int, nnz)
    R = Array(T, nnz)

    if implements_edge_list(g)
        @inbounds @sync @parallel for (i, e)=enumerate(edges(g))
            u, v = e
            ui = vertex_index(u, g)
            vi = vertex_index(v, g)
            I[2i-1] =  i
            J[2i-1] = ui
            R[2i-1] = +1
            I[2i-0] =  i
            J[2i-0] = vi
            R[2i-0] = -1
        end
    elseif implements_adjacency_list(g)
        local idx = 0
        for (ui, u) in enumerate(vertices(g))
            for v in out_neighbors(u, g)
                vi = vertex_index(v, g)
                I[2*idx+1] = idx+1
                J[2*idx+1] = ui
                R[2*idx+1] = +1
                I[2*idx+2] = idx+1
                J[2*idx+2] = vi
                R[2*idx+2] = -1
                idx += 1
            end
        end
    else
        error("g does not implement proper interface.")
    end
    assert (m == idx)
    sparse(I, J, R, m, n)
end
