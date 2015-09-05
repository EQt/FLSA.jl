children = Array{Int}[[], [], [], [1,2,3]]
lb = [1, 2.5, 2, -Inf]
ub = [3.2, 3, 4.5, Inf]
y4 = 2.5

facts("TESTING min_knot! on a simple node (depth 1)") do
    anode = PWLNode(children[4], y4, 4, lb, ub)
    @fact min_knot!(anode) --> 1.0
    @fact min_knot!(anode) --> 2.0
    @fact min_knot!(anode) --> 2.5
    @fact min_knot!(anode) --> 3.0
    @fact min_knot!(anode) --> 3.2
    @fact min_knot!(anode) --> 4.5
    @fact min_knot!(anode) --> Inf
end
