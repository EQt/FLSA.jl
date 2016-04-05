using Convex
using ECOS

function grid_matrix(n1::Integer, n2::Integer)
    n = n1 * n2
    nh = n1*(n2 - 1)    # horizontal
    nv = n2*(n1 - 1)    # vertical
    Dh = spzeros(Int, nh, n)
    for i in 1:n1
        for j in 1:n2-1
            v = (i-1)*(n2-1) + j
            Dh[v, v+0] = -1
            Dh[v, v+1] = +1
        end
    end
    Dv = spzeros(Int, nv, n)
    for i in 1:n1-1
        for j in 1:n2
            v1 = (i-1)*n2 + j
            v2 = (i-0)*n2 + j
            e = v1
            Dv[e, v1] = +1
            Dv[e, v2] = -1
        end
    end
    return vcat(Dh, Dv)
end

function cvx{T<:Real}(y::Array{T,2}; lambda::Real=1.0,
                      solver = ECOSSolver(verbose=false))
    n1, n2 = size(y)
    n = length(y)
    y = reshape(y, n)
    x = Variable(n)
    D = grid_matrix(n1, n2)
    p = Convex.minimize(sum_squares(x - y) + lambda*Convex.norm_1(D*x))
    solve!(p, solver)
    reshape(x.value, (n1,n2))
end
