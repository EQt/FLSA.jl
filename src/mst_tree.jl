max_gap_tree(y::Matrix{Float64}, g::FLSA.ImgGraph; params...) = max_gap_tree(y[:], g; params...)
function max_gap_tree(y::Vector{Float64}, g::FLSA.ImgGraph; c0 = 0.0,
                      max_iter=1,
                      logger = Dict{UTF8String,Any}(),
                      random_tree=false)
    function logg(msg); end
    alpha = c0 * sign(g.D*y[:])
    gap = FLSA.duality_gap(y[:], alpha, g)
    logg("duality gap: $gap")
    dual_obj = FLSA.norm2(y[:] - g.D'*alpha)
    logg("duality obj: $dual_obj")
    for it = 1:max_iter
        if random_tree
            weights = randn(num_edges(g.graph))
        else
            weights = - FLSA.gap_vec(y[:], alpha, g)
        end
        logg("weights")
        mst, wmst = kruskal_minimum_spantree(g.graph, weights)

        logg("tree")
        t = FLSA.subtree(g.graph, mst, 1)
        logg("created subtree")
        z = y[:] - FLSA.non_tree(g.D, mst)'*alpha
        logg("non_tree")
        Lam = fill(Inf, length(y))
        for e in mst
            v, u = source(e), target(e)
            if t.parent[v] == u
                Lam[v] = g.lambda[e.index]
            else
                Lam[u] = g.lambda[e.index]
            end
        end
        logg("lambda subtree")
        x = FLSA.dp_tree(z, Lam , t)
        logg("dp_tree")
        alpha_t = FLSA.dual_tree(z, x, t)
        @debug("gap(tree-part) = $(norm(z - FLSA.tree_part(g.D, mst)' * alpha_t - x))")
        logg("dual_tree: \n$(alpha_t[1:min(5, length(alpha_t))])")
        m = num_edges(g.graph)
        alpha_old = deepcopy(alpha)
        for (i,e) in enumerate(mst)
            alpha[e.index] = alpha_t[i] / g.lambda[e.index]
        end
        # x += z
        gap = FLSA.duality_gap(y[:], alpha, g)
        logg("duality gap: $gap")
        dual_obj = FLSA.norm2(y[:] - g.D'*alpha)
        logg("duality obj: $dual_obj")
    end
end

