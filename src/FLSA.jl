__precompile__()

if VERSION >= v"0.7-"
    using Base.include
    using LinearAlgebra     # for norm
    using Printf
end

module FLSA

const ∞ = Inf
const ℝ = Float64
const Node = Int
const IncMat = SparseMatrixCSC{Float64,Int}

end # module FLSA
