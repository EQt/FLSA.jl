import Base: push!, start, next, done, length, isempty, +

"""Double-ended priority queue: Implementation based on sorted arrays"""
mutable struct DePQ{E}
    elements::Vector{E}
    by::Function
end


if CUSTOM_PRINTER
    import Base.show
    function show{E}(io::IO, pq::DePQ{E})
        Base.showarray(IOContext(io, :limit => true), pq.elements, false)
    end

    import Base.string
    function string{E}(pq::DePQ{E})
        io = IOBuffer()
        show(io, pq)
        takebuf_string(io)
    end
end


@inline front{E}(q::DePQ{E}) = q.elements[1]
@inline back{E}(q::DePQ{E})  = q.elements[end]
@inline pop_front!{E}(q::DePQ{E}) = shift!(q.elements)
@inline pop_back!{E}(q::DePQ{E})  = pop!(q.elements)
@inline function push_front!{E}(q::DePQ{E}, e::E)
    if !isempty(q.elements)
        @assert e.x <= first(q.elements).x "push_front!; e=$e; q=$q"
    end
    unshift!(q.elements, e)
end

@inline function push_back!{E}(q::DePQ{E}, e::E)
    if !isempty(q.elements)
        @assert last(q.elements).x <= e.x "push_back!; e=$e; q=$q"
    end
    push!(q.elements, e)
end


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

@inline length{E}(q::DePQ{E}) = length(q.elements)
@inline isempty{E}(q::DePQ{E}) = isempty(q.elements)

# for convinience (map, enumerate, …)
start{E}(q::DePQ{E}) = start(q.elements)
next{E}(q::DePQ{E}, s) = next(q.elements, s)
done{E}(q::DePQ{E}, s) = done(q.elements, s)
