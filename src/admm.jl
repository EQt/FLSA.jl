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
    size(y,1) == n || error(@sprintf("y has wrong dimension %d (should be %d", size(y,1), n))
    L = D'*D # Laplacian matrix
    @assert size(L, 1) == n
    @assert size(L, 2) == n
    x = copy(y) # initialize x
    b = zeros(m)
    z = zeros(m)
    k::Unsigned = 1 # iteration number
    while k ≤ max_iter
        k += 1
        # update x
        A = eye(L) + μ*L
        c = y + D'*(μ*b - z)
        
        x = conjugate_gradient(A, c, c)
        # update b
        b = soft_threshold(D*x + z/μ, λ/μ, )
        
        # update z
        z += δ * (D*x - b)
        if verbose
            println(norm2(D*x - b))
        end
    end
    x
end

function admm{T<:Number}(y::AbstractMatrix{T},
                         λ::Number = 1.0,
                         δ::Number = 0.1,
                         μ::Number = 0.1;
                         max_iter::Int = 100,
                         verbose::Bool = false)
    n1, n2 = size(y)
    n = n1 * n2
    D = incidence_matrix(y)
    x = admm(reshape(y, n), D, λ, δ, μ; max_iter=max_iter, verbose=verbose)
    return reshape(x, n1, n2)
end
