"""L2 norm squared"""
norm2(x) = (x' * x)[1]

"""The objective function value"""
flsa(x, y, D, λ) = 0.5 * norm2(y-x) + λ*norm(D*x, 1)

"""Input to compute the graph induced fused LASSO signal approximator FLSA"""
type Instance
    y::Vector{Float64}
    lambda::Float64
    graph::AbstractGraph
end

"""Compute the duality gap"""
function duality_gap{T<:Number,I<:Number}(alpha::Vector{T}, lambda::T, y::Vector{T}, D::AbstractMatrix{I})
    psi = D * (D' * alpha - y)
    return lambda * norm(psi, 1) + alpha' * psi
end
