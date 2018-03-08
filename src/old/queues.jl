include("heap.jl")
include("depq.jl")
# include("heap2.jl")

const event_order = Base.Order.By(event_time)
# typealias Q1 SortedMultiDict{Event, Void, typeof(event_order)}
const Q1 = SortedSet{Event,typeof(event_order)}
const Q2 = DePQ{Event}

keys(q::Q1) = q
EventQueue1() = Q1(event_order)
@inline +{E<:Element,O}(q::SortedSet{E,O}, p::SortedSet{E,O}) = merge!(q,p)

EventQueue2() = Q2([], event_time)


macro parent_defines(c)
    local pm = esc(module_parent(current_module()))
    s = string(c)
    return quote
        isdefined($pm, Symbol($s)) &&
            typeof($pm.$c) == Bool &&
            $pm.$c
    end
end


if @parent_defines debug_queues
    import Base: collect, ==
    EventQueue() = (EventQueue1(), EventQueue2())

    ==(e::Event, g::Event) =
        abs(e.x - g.x) <= 1e-11 && e.slope == g.slope && e.offset == g.offset

    pri(v::Vector{Event}) =
        join([@sprintf("%3.2f@%3.2f,%3.2f", e.x, e.slope, e.offset)
              for e in v], "|")
    collect(q::Tuple{Q1,Q2}) = q[1]

    splitq(q) =  collect(keys(q[1])), [e for e in q[2]]
    function print_queues(q)
        q1, q2 = splitq(q)
        "q1=$(pri(q1)), q2=$(pri(q2))"
    end

    function assert_equal(q::Tuple{Q1,Q2})
        q1, q2 = splitq(q)
        @assert length(q1) == length(q2) "q1=$(pri(q1)), q2=$(pri(q2))"
        for i = 1:length(q1)
            @assert q1[i] == q2[i] "q1=$(pri(q1)), q2=$(pri(q2))"
        end
    end
    
    function front(q::Tuple{Q1,Q2})
        e1 = front(q[1])
        e2 = front(q[2])
        @assert e1 == e2
        assert_equal(q)
        e1
    end

    function back(q::Tuple{Q1,Q2})
        e1 = back(q[1])
        e2 = back(q[2])
        @assert e1 == e2
        assert_equal(q)
        e1
    end
    
    function pop_front!(q::Tuple{Q1,Q2})
        e1 = pop_front!(q[1])
        e2 = pop_front!(q[2])
        @assert e1 == e2
        assert_equal(q)
        e1
    end

    function pop_back!(q::Tuple{Q1,Q2})
        e1 = pop_back!(q[1])
        e2 = pop_back!(q[2])
        @assert e1 == e2
        assert_equal(q)
        e1
    end

    function push_front!(q::Tuple{Q1,Q2}, e::Event)
        @debug print_queues(q)
        assert_equal(q)
        push_front!(q[1], e)
        push_front!(q[2], e)
        @debug "push_front!($e)"
        assert_equal(q)
        q
    end

    function push_back!(q::Tuple{Q1,Q2}, e::Event)
        push_back!(q[1], e)
        push_back!(q[2], e)
        assert_equal(q)
        q
    end
    
    function +(q::Tuple{Q1,Q2}, p::Tuple{Q1,Q2})
        r = (q[1]+p[1], q[2]+p[2])
        assert_equal(r)
        r
    end

elseif @parent_defines sortedset
    info("Activating SortedSet")
    include("heap.jl")
    EventQueue = EventQueue1
else
    info("Activating MergedArrays")
    EventQueue = EventQueue2
end
