module TestLine
using FactCheck

debug = true; include("../src/FLSA.jl")

facts("4 line knot") do
    FLSA.@debug("blub")
end

end
