"""L2 norm squared"""
norm2(x) = (x' * x)[1]

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

    r = A*x - b
    p = -r
    r_norm2 = norm2(r)
    k::Unsigned = 1
    while r_norm2 > ɛ && k ≤ max_iter
        α = - (r' * p)[1] / (p' * A * p)[1]
        x += α .* p
        r += α .* A * p
        β = norm2(r)/ r_norm2
        p *=  β
        p -= r
        k += 1
        r_norm2 = norm2(r)
    end
    x
end
