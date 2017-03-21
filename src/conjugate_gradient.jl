using Compat.issymmetric

"""Invoke a fast BLAS routine"""
macro blas(e)
    # println(e)
    if e.head == :.*=
        @assert length(e.args) == 2
        x = e.args[1]
        a = e.args[2]
        ex = :(BLAS.scal!(length($x), $(esc(a)), $(esc(x)), 1))
        # println(ex)
        ex
    else
        # axpy: a*x plus y
        @assert e.head in [:+=, :.+=]
        first = e.args[1]
        rest = e.args[2]
        op = rest.args[1]
        @assert op in [:.*, :*]
        second = rest.args[2]
        third = rest.args[3]
        :(BLAS.axpy!($(esc(second)), $(esc(third)), $(esc(first))))
    end
end


# type CGState
#     r::Vector{Float64}
#     p::Vector{Float64}
#     Ap::Vector{Float64}
# end


"""
conjugate gradient method to solve A*x = b.
Break after max_iter iterations or as soon as ||A*x - b|| < ɛ holds.
"""
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
    k = 1
    Ap = Vector{Float64}(length(r))
    while r_norm2 > ɛ && k ≤ max_iter
        Ap .= A * p
        α = - dot(r, p) / dot(p, Ap)
        @blas x .+= α .* p
        @blas r .+= α .* Ap
        rn = norm2(r)
        @blas p .*= rn / r_norm2
        @blas p += -1.0 * r
        k += 1
        r_norm2 = rn
    end
    x
end
