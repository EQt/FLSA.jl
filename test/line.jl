module TestLine
using FactCheck

debug = true; include("../src/FLSA.jl")

facts("3 line knot") do
    y = [1.0, 2.0, 3.5]
    x = FLSA.dp_line_naive(y, 1.0)
    @fact x --> roughly([2.0, 2.0, 2.5])
end

facts("line: 3 knots, fast") do
    y = [1.0, 2.0, 3.5]
    x = FLSA.dp_line(y, 1.0)
    @fact x --> roughly([2.0, 2.0, 2.5])
end    

end
