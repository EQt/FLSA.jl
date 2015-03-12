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
                                   μ::Number = 0.1)
    m, n = size(D)
    x = copy(y) # initialize x
    b = zeros(n)
    v = zeros(m)
    # TODO: update x
    # TODO: update b
    b = soft_threshold(D*x + v/μ, λ/μ, )
    v += δ * (D*x - b)
end
