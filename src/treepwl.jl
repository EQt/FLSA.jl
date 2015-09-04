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
t = FLSA.subtree(g, mst, (1,1))


type Knot
    islb::Bool
    x::Float64
    i::Int
end

type PWLNode
    events::Vector{Knot}
    a::Int  # index of lowest events, not occured yet
    b::Int  # index of highest events, not occured yet
    slope::Int
    y::Float64
    v::Int  # node index
    lb::Float64
    ub::Float64
    function PWLNode(children, y, v, lb, ub)
        events = [[Knot(true,  lb[c], c) for c in children]
                  [Knot(false, ub[c], c) for c in children]]
        new(events, 1, lenght(events), 1, y, v, -Inf, +Inf)
    end
end

"""Find and extract the next knot from the lower in a node"""
function min_knot!(n::PWLNode, t::TreeSubGraph)
end

"""Find and extract the next knot from the upper in a node"""
function max_knot!(n::PWLNode, t::TreeSubGraph)
end


type PWLTree
    nodes::Vector{PWLNode}
end


"""
function dp_tree_pwl(y::Vector{Float64}, lambda::Float64, t) #t::TreeSubGraph)
    for i in t.dfs_order[end:-1:1]
    end
end
"""

slope = fill(1.0, n)
ub = fill(Inf, n)
lb = fill(-Inf, n)
events = fill(Knot[], n)
lbi = fill(1, n)

for i in t.dfs_order[end:-1:1]
    sort!(events[i], by=k->k.x)
    oldx = y[i]
    lb[i] = y[i] - lambda # prognose
    j = lbi[i]
    while j <= length(events[i])
        j += 1
        knot = events[i][j]
        if lb[i] < knot.x break end
        lb[i] += slope[i]*(knot.x - oldx)
        oldx = knot.x
        if knot.islb
            
        else
            lbi[knot.i] += 1
            knot.i = t.parent[knot.i]
            knot.x = ub[knot.i]
            sort!(events[i], by=k->k.x, lo=j)
            j -= 1
        end
        # TODO: resort events[i]
    end
end

