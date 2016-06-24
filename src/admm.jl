include("conjugate_gradient.jl")

"""Componentwise soft threshold on x by width λ"""
soft_threshold(x, λ) = sign(x) .* max(0, abs(x) - λ)


macro log_admm()
    return quote
        if verbose
            if !haskey(logger, "flsa")
                logger["time"] = []
                logger["flsa"] = []
                logger["ɛ_CG"] = []
            end
            process(x)
            push!(logger["flsa"], flsa(x, y, D))
            push!(logger["time"], time)
            push!(logger["ɛ_CG"], ɛ_CG)
            println(@sprintf("%4d %f %f", k, logger["flsa"][end], norm(x-y, 2)))
        end
    end
end    


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
              max_time::Number = Inf,
              process = x->nothing)
    m, n = size(D)              # Incidence matrix
    L = D'*D                    # Laplacian matrix
    @assert size(y,1) == n @sprintf("size(y)=%d (n=%d)", size(y,1), n)
    @assert size(L) == (n, n)
    x = copy(y)                 # initializcation
    b = zeros(m)
    z = zeros(m)
    k = 1                       # iteration number
    total = 0                   # total runtime
    tic()
    while k ≤ max_iter && total ≤ max_time
        total += toq()                          # stop cloc for logging
        @log_admm
        tic()
        A = eye(L) + μ*L                        # lhs matrix
        c = y + D'*(μ*b - z)                    # rhs
        x = conjugate_gradient(A, c, x; ɛ=ɛ_CG) # update x
        ɛ_CG *= ɛ_c                             # update accuracy
        b = soft_threshold(D*x + z/μ, 1/μ)      # update b
        z += δ * (D*x - b)                      # update z
        k += 1
    end
    x
end

"""For convenience…"""
admm(y::Matrix{Float64}, g::ImgGraph; params...) =
    reshape(admm(y[:], g.D; params...), g.n1, g.n2)
