# using Devectorize

macro logitfista()
    return quote
        if verbose
            if !haskey(logger, "flsa")
                logger["flsa"] = []
                logger["time"] = []
                logger["gap"]  = []
            end
            x = y - D'*α
            push!(logger["flsa"], flsa(x, y, D))
            push!(logger["time"], time)
            push!(logger["gap"], -1.0)
                  # duality_gap(α, y, D))
            println(@sprintf("%4d %f %f", k,
                             logger["flsa"][end], logger["gap"][end]))
        end
    end
end



"""
Compute FLSA by FAST GRADIENT PROJECTION.
Also called Fast Iterative Shrinkage/Thresholding Algorithm.

Notice that calculation is done in normed dual variables.

`D` : *weighted* oriented incidence matrix
"""
function fista(y::Vector{Float64},
               D::IncMat;
               L::Float64 = 8,
               max_iter::Int = 100,
               verbose::Bool = false,
               logger = Dict{UTF8String,Any}(),
               max_time::Float64 = Inf)
    m, n = size(D)
    size(y,1) == n ||
      error(@sprintf("y has wrong dimension %d (should be %d", size(y,1), n))

    prox(x) = clamp(x, -1.0, +1.0)
    grad(α) = D*(D'*α - y)              # gradient
    pL(α) = prox(α - 1/L*grad(α))

    tic()
    total = 0
    α = β = sign(D * y)
    t = 1
    k = 1
    while k <= max_iter+1 && total ≤ max_time
        total += (time = toq())
        @logitfista
        if k == max_iter break end
        tic()
        α₀ = α
        α = pL(β)
        t₁ = (1 + sqrt(1 + 4t^2))/2
        β = α + (t - 1)/t₁ * (α - α₀)
        t = t₁
        k += 1
    end
    return y - D'*α
end

"""For convenience…"""
function fista{T<:Number}(y::AbstractMatrix{T},
                                    D::SparseMatrixCSC{Float64,Int};
                                    params...)
    n1, n2 = size(y)
    x = fista(reshape(y, n1*n2), D; params...)
    return reshape(x, n1, n2)
end
