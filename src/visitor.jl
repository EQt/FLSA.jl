"""Contains all information for the calculation"""
type DPVisitor <: AbstractGraphVisitor
    y::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    df::Vector{PWL}
    function DPVisitor(y::Vector{Float64})
        n = length(y)
        lb = fill(-Inf, n)
        ub = fill(+Inf, n)
        df = [PWL([]) for i in y]
        new(y, lb, ub, df)
    end
end

"""For testing the graph traversal mechanism"""
function Graphs.close_vertex!(vis::DPVisitor, v)
    println("processing $v")
end
