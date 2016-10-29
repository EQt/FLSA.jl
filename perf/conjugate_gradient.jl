norm2(x::Vector{Float64}) = dot(x,x)

square{MatT<:AbstractMatrix{Float64}}(A::MatT, p::Vector{Float64}) =
    (p' * A * p)[1]

# include("../src/conjugate_gradient.jl")
function conjugate_gradient{MatT<:AbstractMatrix{Float64}}(A::MatT,
                                                           b::Vector{Float64},
                                                           x::Vector{Float64};
                                                           ɛ::Float64 = 1e-6,
                                                           max_iter::Int = 100)
    size(A,1) == size(A,2) || error("Matrix A has to be a square matrix")
    issymmetric(A) || error("Matrix A has to be a symmetric matrix")

    r = A*x - b
    p = -r
    ɛ = max(ɛ*norm(A, Inf) * norm(b, Inf), 1e-6)
    r_norm2 = norm2(r)
    k::Unsigned = 1
    while r_norm2 > ɛ && k ≤ max_iter
        α = - dot(r, p) / square(A, p)
        x .+= α .* p
        r .+= α .* A * p
        β = norm2(r)/ r_norm2
        p .*=  β
        BLAS.axpy!(-1, r, p)
        k += 1
        r_norm2 = norm2(r)
    end
    x
end


n = 100
A = Matrix(Symmetric(rand(n,n)))
@assert issymmetric(A)
b = rand(n)
x = zeros(b)

conjugate_gradient(A, b, x)

