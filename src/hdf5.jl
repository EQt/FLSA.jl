using HDF5

# file_name = Pkg.dir("FLSA", "examples", "example.h5")
function read_h5(file_name::AbstractString)
    h5open(file_name, "r") do io
        y  = read(io, "nodes/input")
        mu = read(io, "nodes/weight")
        I  = read(io, "edges/head")
        J  = read(io, "edges/tail")
        V  = read(io, "edges/weight")
        n  = length(y)
        m  = length(I)
        @assert m == length(J)
        @assert m == length(V)
        D = sparse([collect(1:m); collect(1:m)], Int[I;J;], [V; -V], m, n)
        return mu, D
    end
end
