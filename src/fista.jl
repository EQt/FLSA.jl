"""Clip x within [l,u]"""
clip(x, l, u) = min(u, max(l, x))

"""
Compute FLSA by FAST ITERATIVE SHRINKAGE/THRESHOLDING ALGORITHM.
"""
function fista{T<:Number,I<:Number}(y::Vector{T},
                                    D::AbstractMatrix{I},
                                    λ::Number = 1.0;
                                    L = 8,
                                    max_iter::Int = 100,
                                    verbose::Bool = false,
                                    logger = Dict{String,Any}())
    m, n = size(D)
    size(y,1) == n ||
      error(@sprintf("y has wrong dimension %d (should be %d", size(y,1), n))

    prox(x) = clip(x, -λ, +λ)
    grad(α) = D*(D'*α - y)              # gradient
    pL(α) = prox(α - 1/L*grad(α))

    α = β = λ * sign(D * y)
    t = 1
    k = 1
    while k <= max_iter
        tic()
        α₀ = α
        α = pL(β)
        t₁ = (1 + sqrt(1 + 4t^2))/2
        β = α + (t - 1)/t₁ * (α - α₀)
        t = t₁
        time = toq()
        if verbose
            if !haskey(logger, "flsa")
                logger["flsa"] = {}
                logger["time"] = {}
            end
            x = y - D'*α
            push!(logger["flsa"], flsa(x, y, D, λ))
            push!(logger["time"], time)
            println(@sprintf("%4d %f", k, logger["flsa"][end]))
        end
        k += 1
    end
    return y - D'*α
end

"""For convenience…"""
function fista{T<:Number,I<:Number}(y::AbstractMatrix{T},
                                    D::AbstractMatrix{I},
                                    λ::Number = 1.0;params...)
    n1, n2 = size(y)
    x = fista(reshape(y, n1*n2), D, λ; params...)
    return reshape(x, n1, n2)
end

