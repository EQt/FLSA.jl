module TestLine
using FactCheck
using FLSA

reload("FLSA")

facts("3 line knot") do
    y = [1.0; 2.0; 3.5]
    x = FLSA.dp_line_naive(y, 1.0)
    @fact x --> roughly([2.0; 2.0; 2.5])
end

facts("line: 3 knots, fast") do
    y = [1.0, 2.0, 3.5]
    x = FLSA.dp_line(y, 1.0)
    @fact x --> roughly([2.0; 2.0; 2.5])
end

facts("line: 5 knots same") do
    y = [1.0; 4.0; 3.5; 3.8; 4.1]
    z = FLSA.dp_line_naive(y, 1.0)
    x = FLSA.dp_line(y, 1.0)
    @fact x --> roughly(z)
end

facts("line: 3 rand knots same") do
    srand(13)
    y = rand(20)
    z = FLSA.dp_line_naive(y, 0.3)
    x = FLSA.dp_line(y, 0.3)
    d = 3
    @fact round(x,d) --> roughly(round(z,d))
end

end
