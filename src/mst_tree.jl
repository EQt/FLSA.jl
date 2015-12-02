max_gap_tree(y::Matrix{Float64}, g::FLSA.ImgGraph; params...) = max_gap_tree(y[:], g; params...)
function max_gap_tree(y::Vector{Float64}, g::FLSA.ImgGraph; c0 = 0.0,
                      max_iter=1,
                      logger = Dict{UTF8String,Any}())
    global alpha = c0 * sign(grid.D*y[:])
    # global alpha = zeros(num_edges(grid.graph))
    gap = FLSA.duality_gap(y[:], alpha, grid)
    logg("duality gap: $gap")
    dual_obj = FLSA.norm2(y[:] - grid.D'*alpha)
    logg("duality obj: $dual_obj")
    # imshow(y; title="noisy")
    for it = 1:max_iter
        global weights = - FLSA.gap_vec(y[:], alpha, grid)
        # global weights = randn(num_edges(grid.graph))
        logg("computed weights")
        mst, wmst = kruskal_minimum_spantree(grid.graph, weights)
        global mst, wmst
        logg("computed tree")
        global t = FLSA.subtree(grid.graph, mst, 1)
        logg("created subtree")
        global z = y[:] - FLSA.non_tree(grid.D, mst)'*alpha
        logg("non_tree")
        # imshow(reshape(z, n1, n2))
        global Lam = fill(Inf, length(y))
        for e in mst
            v, u = source(e), target(e)
            if t.parent[v] == u
                Lam[v] = grid.lambda[e.index]
            else
                Lam[u] = grid.lambda[e.index]
            end
        end
        logg("lambda subtree")
        global x = FLSA.dp_tree(z, Lam , t)
        logg("dp_tree")
        global alpha_t = FLSA.dual_tree(z, x, t)
        warn(norm(z - FLSA.tree_part(grid.D, mst)' * alpha_t - x))
        logg("dual_tree: \n$(alpha_t[1:min(5, length(alpha_t))])")
        global m = num_edges(grid.graph)
        global alpha_old = deepcopy(alpha)
        for (i,e) in enumerate(mst)
            alpha[e.index] = alpha_t[i] / grid.lambda[e.index]
        end
        # x += z
        gap = FLSA.duality_gap(y[:], alpha, grid)
        logg("duality gap: $gap")
        dual_obj = FLSA.norm2(y[:] - grid.D'*alpha)
        logg("duality obj: $dual_obj")
        imshow(reshape(x, n1, n2); title="tree $it")
    end
end

