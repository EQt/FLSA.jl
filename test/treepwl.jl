children = Array{Int}[[], [], [], [1,2,3]]
lb = [1, 2.5, 2, -Inf]
ub = [3.2, 3, 4.5, Inf]
y4 = 2.5
anode = PWLNode(children[4], y4, 4, lb, ub)

@assert (x = min_knot!(anode)) == 1.0 "x = $x"
@assert (x = min_knot!(anode)) == 2.0 "x = $x"
