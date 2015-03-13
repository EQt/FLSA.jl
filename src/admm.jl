include("conjugate_gradient.jl")

"""Componentwise soft threshold on x by width λ"""
soft_threshold(x, λ) = sign(x) .* max(0, abs(x) - λ)

"""
Solve the FLSA by the ALTERNATING DIRECTION METHOD OF MULTIPLIERS.
"""
function admm{T<:Number,I<:Number}(y::Vector{T},
                                   D::AbstractMatrix{I},
                                   λ::Number = 1.0,
                                   δ::Number = 0.1,
                                   μ::Number = 0.1;
                                   max_iter::Int = 100,
                                   verbose::Bool = false)
    m, n = size(D)
    L = D'*D # Laplacian matrix
    x = copy(y) # initialize x
    b = zeros(n)
    v = zeros(m)
    k::Unsigned = 1 # iteration number
    while k ≤ max_iter
        k += 1
        # update x
        A = eye(L) + μ*L
        c = y + D'*(μ*b - z)
        
        x = conjugate_gradient(A, c, c)
        # update b
        b = soft_threshold(D*x + v/μ, λ/μ, )
        
        # update v
        v += δ * (D*x - b)
    end
end
