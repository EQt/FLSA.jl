"""More commonly used name"""
clip(x::Float64, lo::Float64, hi::Float64) = Base.clamp(x, lo, hi)
clip(x::Vector{Float64}, lo::Float64, hi::Float64) = Base.clamp(x, lo, hi)


"""L2 norm squared"""
norm2(x) = dot(x, x)

"""The objective function value"""
flsa(x, y, D, λ=1.0) = 0.5 * norm2(y-x) + λ*norm(D*x, 1)


"""Input to compute the graph induced fused LASSO signal approximator FLSA"""
type Instance
    y::Vector{Float64}
    lambda::Float64
    graph::AbstractGraph
end

"""Compute the duality gap"""
function duality_gap{T<:Number,I<:Number}(alpha::Vector{T}, lambda::T, y::Vector{T}, D::AbstractMatrix{I})
    psi = D * (D' * alpha - y)
    return lambda * norm(psi, 1) + dot(alpha, psi)
end

function gap_vec(y, alpha, D::IncMat)
    @assert minimum(alpha) >= -1.0 - 1e-9 @val(minimum(alpha))
    @assert maximum(alpha) <= +1.0 + 1e-9 @val(maximum(alpha))
    x = y - D' * alpha
    g = - D*x
    return (alpha .* g) + abs(g)
end
