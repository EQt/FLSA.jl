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
                                   verbose::Bool = false,
                                   logger = Dict{String, Any}())
    m, n = size(D)
    size(y,1) == n ||
       error(@sprintf("y has wrong dimension %d (should be %d", size(y,1), n))

    L = D'*D # Laplacian matrix
    @assert size(L, 1) == n
    @assert size(L, 2) == n
    x = copy(y) # initialize x
    b = zeros(m)
    z = zeros(m)
    k= 1 # iteration number
    tic()
    while k ≤ max_iter
        time = toq()
        if verbose
            if !haskey(logger, "flsa")
                logger["time"] = {}
                logger["flsa"] = {}
            end
            push!(logger["flsa"], flsa(x, y, D, λ))
            push!(logger["time"], time)
            println(@sprintf("%4d %f", k, logger["flsa"][end]))
        end
        tic()
        # update x
        A = eye(L) + μ*L
        c = y + D'*(μ*b - z)
        
        x = conjugate_gradient(A, c, x)
        # update b
        b = soft_threshold(D*x + z/μ, λ/μ, )
        
        # update z
        z += δ * (D*x - b)
        k += 1
    end
    x
end

"""For convenience…"""
function admm{T<:Number}(y::AbstractMatrix{T},
                         D::AbstractMatrix{T},
                         λ::Number = 1.0,
                         δ::Number = 0.1,
                         μ::Number = 0.1; params...)
    n1, n2 = size(y)
    x = admm(reshape(y, n1*n2), D, λ, δ, μ; params...)
    return reshape(x, n1, n2)
end


"""For convenience (see first admm function)"""
admm{T<:Number}(y::AbstractMatrix{T}, λ::Number = 1.0, δ::Number = 0.1, μ::Number = 0.1; params...) =
    admm(y, grid_graph(size(y)...), λ, δ, μ; params...)
