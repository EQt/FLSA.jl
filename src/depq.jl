import Base: push!, start, next, done, length, isempty


"""Double-ended priority queue: Implementation based on sorted arrays"""
type DePQ{E}
    elements::Vector{E}
    by::Function
end

typealias MDEPQ{E} DePQ{E}

EventQueue() = MDEPQ{Event}([], e->e.x)

front{E}(q::DePQ{E}) = q.elements[1]
back{E}(q::DePQ{E})  = q.elements[end]
pop_front!{E}(q::DePQ{E}) = shift!(q.elements)
pop_back!{E}(q::DePQ{E})  = pop!(q.elements)
push_front!{E}(q::DePQ{E}, e::E) = unshift!(q.elements, e)
push_back!{E}(q::DePQ{E}, e::E) = push!(q.elements, e)

+{E}(q::DePQ{E}, p::DePQ{E}) = merge!(p, q)
function merge!{E}(q::DePQ{E}, p::DePQ{E})
    append!(q.elements, p.elements)
    sort!(q.elements, by=q.by)
    return q
end

function push!{E}(q::DePQ{E}, e::E)
    push!(q.elements, e)
    sort!(q.elements, by=q.by)
end

length{E}(q::DePQ{E}) = length(q.elements)
isempty{E}(q::DePQ{E}) = length(q) == 0

# for convinience (map, enumerate, …)
start{E}(q::DePQ{E}) = start(q.elements)
next{E}(q::DePQ{E}, s) = next(q.elements, s)
done{E}(q::DePQ{E}, s) = done(q.elements, s)
