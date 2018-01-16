if isfile(Pkg.dir("HDF5", "src", "HDF5.jl"))

using HDF5

# file_name = Pkg.dir("FLSA", "examples", "example.h5")
function read_h5(file_name::AbstractString)
    h5open(file_name, "r") do io
        y  = read(io, "nodes/input")
        mu = read(io, "nodes/weight")
        H  = read(io, "edges/head")
        T  = read(io, "edges/tail")
        L  = read(io, "edges/weight")
        n  = length(y)
        m  = length(H)
        @assert m == length(T)
        @assert m == length(L)
        @assert n == length(mu)
        @assert 1 <= minimum(T) "minimum(edges/tail) = $(minimum(T))"
        @assert 1 <= minimum(H) "minimum(edges/head) = $(minimum(H))"
        @assert maximum(T) <= n "n=$(n), maximum(edges/tail)=$(maximum(T))"
        @assert maximum(H) <= n "n=$(n), maximum(edges/head)=$(maximum(H))"
        I = [collect(1:m); collect(1:m)]
        J = Int[H;T]
        V = [L; -L]
        D = sparse(I, J, V, m, n)
        return y, mu, D
    end
end

end
