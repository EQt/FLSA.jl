"Compute the derivate for a line and exentds bounds accordingly"
function dp_line(y::Vector{Number}, bounds)
    error("Not implemented yet")
end

type DPVisitor <: AbstractGraphVisitor
    y::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    df::Vector{PWL}
    function DPVisitor(y::Vector{Float64})
        n = length(y)
        lb = fill(-Inf, n)
        ub = fill(+Inf, n)
        df = [PWL([Knot(0,yi)]) for yi in y]
        println("Constructing visitor")
        new(y, lb, ub, df)
    end
end

function Graphs.close_vertex!(vis::DPVisitor, v)
    println("processing $v")
end


