using Devectorize

macro logitfista()
    return quote
        if verbose
            if !haskey(logger, "flsa")
                logger["flsa"] = {}
                logger["time"] = {}
                logger["gap"] = {}
            end
            x = y - D'*α
            push!(logger["flsa"], flsa(x, y, D, λ))
            push!(logger["time"], time)
            push!(logger["gap"], duality_gap(α, λ, y, D))
            println(@sprintf("%4d %f %f", k,
                             logger["flsa"][end], logger["gap"][end]))
        end
    end
end


"""
Compute FLSA by FAST ITERATIVE SHRINKAGE/THRESHOLDING ALGORITHM.
"""
function fista{MT<:AbstractMatrix}(y::Vector{Float64},
                                   D::MT,
                                   λ::Float64 = 1.0;
                                   L::Float64 = 8,
                                   max_iter::Int = 100,
                                   verbose::Bool = false,
                                   logger = Dict{String,Any}(),
                                   max_time::Float64 = Inf)
    m, n = size(D)
    size(y,1) == n ||
      error(@sprintf("y has wrong dimension %d (should be %d", size(y,1), n))

    prox(x) = clamp(x, -λ, +λ)
    grad(α) = D*(D'*α - y)              # gradient
    pL(α) = prox(α - 1/L*grad(α))

    tic()
    total = 0
    α = β = λ * sign(D * y)
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
function fista{T<:Number,I<:Number}(y::AbstractMatrix{T},
                                    D::AbstractMatrix{I},
                                    λ::Number = 1.0;params...)
    n1, n2 = size(y)
    x = fista(reshape(y, n1*n2), D, λ; params...)
    return reshape(x, n1, n2)
end

