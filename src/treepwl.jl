"""
Faster implementation
~~~~~~~~~~~~~~~~~~~~~
"""

"""Record what is happening, when a knot of the PWL is hit"""
type Event
    s::Int          # source node
    t::Int          # target node
    x::Float64      # position
    offset::Float64 # delta offset
    slope::Float64  # delta slope
end

forecast(e::Event, c) = (c - e.offset)/e.slope
set_forecast(e::Event, c) = (e.x = forecast(e, c))


"""Manage the events of a node"""
type PWLNode
    lb::Float64             # lower bound (computed by create_min_event)
    ub::Float64             # upper bound (computed by create_max_event)
    pq::PDeQue{Event}       # events
    PWLNode() = new(-Inf, +Inf, PDeQue{Event}([], e->e.x))
end


type PWLTree
    nodes::Vector{PWLNode}
    children::Vector{Vector{Int}}
    pre_order::Vector{Int}
    parent::Vector{Int}
    root::Int
    y::Vector{Float64}
    lam::Function

    function PWLTree(parents::Vector{Int}, root::Int, y::Vector{Float64},
                     lambda=i->1.0)
        n = length(parents)
        children = [Int[] for i=1:n]
        for (v,p) in enumerate(parents)
            if v != root push!(children[p], v) end
        end
        nodes = [PWLNode() for i in 1:n]
        pre_order = zeros(Int, n)
        stack = [root]
        nr = 1
        while !isempty(stack)
            v = pop!(stack)
            pre_order[nr] = v
            nr += 1
            append!(stack, children[v])
        end
        return new(nodes, children, pre_order, parents, root, y, lambda)
    end
    PWLTree(t::ITreeSubGraph, y::Vector{Float64}, lambda=i->1.0) =
        PWLTree(t.parent, t.root, y, lambda)
end


find_min(n) = front(n.events)
find_max(n) = back(n.events)

function step_min(t, ek)
    @debug "step_min($(ek.t)): $ek"
    n = t.nodes[ek.t]
    ekk = pop_front!(n)
    @debug "step_min($(ek.t)): ekk = $ekk (will be deleted)"
    @debug "step_min(): Going from $(ek.t) to $(ekk.t)"
    ek.t = ekk.t
    ek.slope  += ekk.slope
    ek.offset += ekk.offset
    ek.x = ekk.x
    pq = t.nodes[ek.s].events
    ekk = pop_front!(pq)
    if abs(ek.slope) <= 1e-6
        @debug "step_min(): trying to delete ek = $ek"
        @assert ek == ekk
        return
    else
        push!(pq, ek)
    end
end


function lower_event!(t, i)
end


function upper_event!(t, i)
end


function print_min(t)
end


function print_max(t)
end


function print_tree(t)
end


function forward_dp_treepwl(t)
    for i in t.pre_order[end:-1:1]
        lower_event!(t, i)
        upper_event!(t, i)
        print_tree(t)
        print_min(t)
        print_max(t)
    end
end
