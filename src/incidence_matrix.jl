"""Compute the oriented incidence matrix of g and store in a sparse matrix of type T"""
function incidence_matrix{T<:Number}(g::AbstractGraph,  ::Type{T} = Int)
    is_directed(g) && error("g must be undirected.")

    local n = num_vertices(g)
    local m = num_edges(g)
    local nnz = 2*m
    I = Array(Int, nnz)
    J = Array(Int, nnz)
    R = Array(T, nnz)

    if implements_edge_list(g)
        if isa(g, GridGraph)
            local n1 = num_rows(g)
            vertex_index(v, g) = n1 * (v[2]-1) + v[1]
        end
        # @inbounds @sync @parallel for e=edges(g)
        @inbounds for e=edges(g)
            u, v, i = source(e), target(e), e.index
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
        assert(m == idx)
    else
        error("g does not implement proper interface.")
    end
    sparse(I, J, R, m, n)
end

incidence_matrix{T}(n1::Int, n2::Int, ::Type{T} = Int) =
    incidence_matrix(grid_graph(n1, n2), T)


incidence_matrix{T<:Number,I<:Number}(y::AbstractMatrix{T}, ::Type{I} = Int) =
    incidence_matrix(size(y)..., I)

