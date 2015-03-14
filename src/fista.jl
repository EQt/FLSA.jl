"""
Compute FLSA by FAST ITERATIVE SHRINKAGE/THRESHOLDING ALGORITHM.
"""
function fista{T<:Number,I<:Number}(y::Vector{T},
                                    D::AbstractMatrix{I},
                                    λ::Number = 1.0;
                                    L = 8,
                                    max_iter::Int = 100,
                                    verbose::Bool = false)
    prox1(x) = soft_threshold(x, 1/L)   # prox operator
    grad(α) = D*(D'*α - y)              # gradient
    pL(α) = prox1(α - 1/L*grad(α))

    m, n = size(D)
    size(y,1) == n || error(@sprintf("y has wrong dimension %d (should be %d", size(y,1), n))
    α = β = λ * sign(D * y)
    t = 1
    k = 1
    while k <= max_iter
        α₀ = α
        α = pL(β)
        t₁ = (1 + sqrt(1 + 4t^2))/2
        β = α + (t - 1)/t₁ * (α - α₀)
        t = t₁
        k += 1
        if verbose
            info(@sprintf("t = %f\t||α - α₀|| = %f", t, norm2(α - α₀)))
        end
    end
    return y - D'α 
end
