"""Piecewise linear function over a tree"""
using FLSA
using Graphs

srand(42)

lambda = 0.1
n1, n2 = 4, 2
g = grid_graph(n1, n2)
v = vertices(g)
n = length(v)
w = rand(size(g.edges))
y = rand(n1*n2)

mst, wmst = kruskal_minimum_spantree(g, w)
tm = FLSA.subtree(g, mst, (1,1))


type Knot
    islb::Bool
    x::Float64
    i::Int
end

slope = fill(1.0, n)
ub = fill(Inf, n)
lb = fill(-Inf, n)
events = fill(Knot[], n)

for i in t.dfs_order[end:-1:1]
    events[i] = sort!([[Knot(true,  lb[c], c) for c in t.children[i]]
                       [Knot(false, ub[c], c) for c in t.children[i]]])
    x = y[i] - lambda

end

function dp_tree_pwl(y::Vector{Float64}, lambda::Float64, t::TreeSubGraph)
    for i in t.dfs_order[end:-1:1]
        # compute lb[i]
        lbis = deepcopy(t.children[i])
        sort!(lbis, by=c->lb[c])
        if lb[first(lbis)] > y[i] - lambda
            lb[i] = y[i] - lambda
        else
            lbx[i] = lb[first(lbis)]
        end
    end
end
