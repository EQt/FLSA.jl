using Compat.issymmetric

"""
conjugate gradient method to solve A*x = b.
Break after max_iter iterations or as soon as ||A*x - b|| < ɛ holds.
"""
function conjugate_gradient{T<:Number, R<:Real}(A::AbstractMatrix{T},
                                                b::Vector{R},
                                                x::Vector{R};
                                                ɛ::Real = 1e-6,
                                                max_iter::Integer = 100)
    size(A,1) == size(A,2) || error("Matrix A has to be a square matrix")
    issymmetric(A) || error("Matrix A has to be a symmetric matrix")

    r = A*x - b
    p = -r
    ɛ = max(ɛ*norm(A, Inf) * norm(b, Inf), 1e-6)
    r_norm2 = norm2(r)
    k::Unsigned = 1
    while r_norm2 > ɛ && k ≤ max_iter
        α = - (r' * p)[1] / (p' * A * p)[1]
        @inplace x += α .* p
        @inplace r += α .* A * p
        β = norm2(r)/ r_norm2
        p *=  β
        BLAS.axpy!(-1, r, p)
        k += 1
        r_norm2 = norm2(r)
    end
    x
end
