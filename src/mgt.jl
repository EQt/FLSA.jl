max_gap_tree(y::Matrix{Float64}, g::FLSA.ImgGraph; params...) =
    reshape(max_gap_tree(y[:], g; params...), g.n1, g.n2)

ONE_FUNCTION = i -> 1.0

function max_gap_tree(y::Vector{Float64}, g::FLSA.ImgGraph;
                      c0::Real = 0.0,
                      mu_f::Function = ONE_FUNCTION,
                      mu::Vector{Float64} = Vector{Float64}(),
                      alpha::Vector{Float64} = Vector{Float64}(),
                      max_iter::Integer=1,
                      logger::LoggerT = LoggerT(),
                      random_tree::Bool=false,
                      abs_tree::Bool=false,
                      verbose::Bool=true,
                      process::Function=x->nothing,
                      dprocess::Function=α->nothing,
                      assert_decreasing::Bool=false)
    function logg(msg...); end
    x = y
    if length(alpha) <= 0
        alpha = c0 * sign(g.D*y[:])
        logg("sign alpha")
    else
        x = y[:] - g.D' * alpha
        logg("predfined x")
    end
    total = 0.0
    tic()
    for it = 0:max_iter
        if verbose
            time = toq()
            total += time
            _field(logger, "time", time)
            _field(logger, "flsa", flsa(x, y, g.D))
            _field(logger, "gap", FLSA.duality_gap(y[:], alpha, g))
            _field(logger, "dual", FLSA.dual_obj(alpha, y, g.D))
            process(x)
            dprocess(alpha)
            println(@sprintf("%4d %f %f %f",
                             it,
                             logger["flsa"][end],
                             logger["dual"][end],
                             logger["gap"][end]))
            if assert_decreasing && length(logger["flsa"]) >= 2
                @assert logger["flsa"][end] <= logger["flsa"][end-1]
            end
        end
        it >= max_iter && break

        tic()
        if random_tree
            weights = randn(num_edges(g.graph))
        elseif abs_tree
            weights = abs(alpha)
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
        logg("sub_lambda")

        x = if length(mu) > 0
            FLSA.dp_tree(z, Lam, mu, t)
        else
            FLSA.dp_tree(z, i->Lam[i], mu_f, t)
        end
        logg("dp_tree")
        alpha_t = FLSA.dual_tree(z, x, t)
        @debug("gap(tree-part) = $(norm(z - FLSA.tree_part(g.D, mst)' * alpha_t - x))")
        logg("dual_tree: \n$(alpha_t[1:min(5, length(alpha_t))])")

        for (i,e) in enumerate(mst)
            alpha[e.index] = alpha_t[i] / g.lambda[e.index]
        end
    end
    return x
end

