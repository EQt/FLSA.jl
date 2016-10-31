import Base.norm

"""More commonly used name"""
clip(x::Float64, lo::Float64, hi::Float64) = Base.clamp(x, lo, hi)
clip(x::Vector{Float64}, lo::Float64, hi::Float64) = Base.clamp(x, lo, hi)


"""L2 norm squared"""
norm2(x::Vector{Float64}) = dot(x, x)
norm2(x::Matrix{Float64}) = norm2(x[:])


"""The objective function value"""
flsa(x::Matrix{Float64}, y::Matrix{Float64}, D::IncMat) =
    flsa(vec(x), vec(y), D)
flsa(x::Vector{Float64}, y::Vector{Float64}, D::IncMat) =
    0.5 * norm2(y-x) + norm(D*x, 1)
flsa0(x::Vector{Float64}, y::Vector{Float64}, D::IncMat, mu::Vector{Float64}) =
    0.5dot(mu, (y-x).^2) + norm(D*x, 1)


"""Input to compute the graph induced fused LASSO signal approximator FLSA"""
type Instance
    y::Vector{Float64}
    lambda::Float64
    graph::AbstractGraph
end

"""Compute the duality gap"""
function duality_gap{T<:Number,I<:Number}(alpha::Vector{T},
                                          lambda::T, y::Vector{T},
                                          D::AbstractMatrix{I})
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


"""Dual objective function"""
dual_obj(alpha::Vector{Float64}, y::Vector{Float64}, D::IncMat) =
    norm2(y - D'*alpha)

dual_obj(alpha::Vector{Float64}, y::Matrix{Float64}, D::IncMat) =
    dual_obj(alpha, vec(y), D)


"""Compute normed histograms"""
normed1(x::Vector) = x /= sum(x)
"""Compute normed histogram"""
nhist(x::Vector) = normed1(hist(x)[2])
nhist(x, e) = normed1(hist(x, e)[2])
export nhist


typealias LoggerT Dict{String,Vector{Float64}}

"""Logg iteration information"""
function _field(logger::LoggerT, name::String, value::Float64)
    if !haskey(logger, name)
        logger[name] = Float64[]
    end
    push!(logger[name], value)
end
