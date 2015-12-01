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


function gap_vec(y, alpha, grid)
    a = alpha ./ grid.lambda
    @assert minimum(alpha) >= -1.0
    @assert maximum(alpha) <= +1.0
    x = y - grid.D' * a
    g = - grid.D*x
    return (a .* g) + abs(g)
end

duality_gap(y, alpha, grid) = sum(gap_vec(y, alpha, grid))
