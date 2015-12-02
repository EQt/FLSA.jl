include("conjugate_gradient.jl")

"""Componentwise soft threshold on x by width λ"""
soft_threshold(x, λ) = sign(x) .* max(0, abs(x) - λ)

"""
Solve the FLSA by the ALTERNATING DIRECTION METHOD OF MULTIPLIERS.
"""
function admm(y::Vector{Float64},
              D::IncMat;
              μ::Number = 0.5,
              δ::Number = 0.5,
              ɛ_CG::Real= 0.1,
              ɛ_c::Real = 0.5,
              max_iter::Int = 100,
              verbose::Bool = false,
              logger = Dict{UTF8String, Any}(),
              max_time::Number = Inf)
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
    total = 0
    tic()
    while k ≤ max_iter && total ≤ max_time
        time = toq()
        total += time
        if verbose
            if !haskey(logger, "flsa")
                logger["time"] = []
                logger["flsa"] = []
                logger["ɛ_CG"] = []
            end
            push!(logger["flsa"], flsa(x, y, D))
            push!(logger["time"], time)
            push!(logger["ɛ_CG"], ɛ_CG)
            println(@sprintf("%4d %f", k, logger["flsa"][end]))
        end
        tic()
        # update x
        A = eye(L) + μ*L
        c = y + D'*(μ*b - z)
        x = conjugate_gradient(A, c, x; ɛ=ɛ_CG)
        ɛ_CG *= ɛ_c
        # update b
        b = soft_threshold(D*x + z/μ, λ/μ)
        
        # update z
        z += δ * (D*x - b)
        k += 1
    end
    x
end

"""For convenience…"""
admm(y::Matrix{Float64}, g::ImgGraph; params...) =
    reshape(admm(y[:], g.D; params...), g.n1, g.n2)
