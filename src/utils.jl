import Base.norm

"""More commonly used name"""
clip(x::Float64, lo::Float64, hi::Float64) = Base.clamp(x, lo, hi)
clip(x::Vector{Float64}, lo::Float64, hi::Float64) = Base.clamp(x, lo, hi)


"""L2 norm squared"""
norm2(x::Vector{Float64}) = dot(x, x)

"""The objective function value"""
flsa(x, y::Vector{Float64}, D::IncMat) =
    0.5 * norm2(y-x) + norm(D*x, 1)


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

function gap_vec(y::Vector{ℝ}, alpha::Vector{ℝ}, D::IncMat)
    m, n = size(D)
    @assert n == length(y) y, n, size(y)
    @assert m == length(alpha) m, size(alpha)
    @assert minimum(alpha) >= -1.0 - 1e-9 @val(minimum(alpha))
    @assert maximum(alpha) <= +1.0 + 1e-9 @val(maximum(alpha))
    x = y - D' * alpha
    g = - D*x
    return (alpha .* g) + abs(g)
end
# convinience
gap_vec(y::Matrix{ℝ}, alpha::Vector{ℝ}, D::IncMat) =
    gap_vec(reshape(y, prod(size(y))), alpha, D)

duality_gap(y::Vector{ℝ}, alp::Vector{ℝ}, D::IncMat) = sum(gap_vec(y, alp, D))

"""Overload the `+=` operator for arrays"""
macro inplace(ex)
    if ex.head == :+=
        # broadcast(+, $(ex.args[1]), $(ex.args[2]), $(ex.args[2]))
        :(BLAS.axpy!(1.0, $(ex.args[2]), $(ex.args[1])))
    else
        ex
    end
end
