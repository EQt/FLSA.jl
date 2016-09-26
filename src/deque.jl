import DataStructures: back, front, Deque
import Base: +, merge!

typealias DeQue Deque


"""More convinient way of denoting the functions"""
pop_front!{T}(q::DeQue{T}) = shift!(q)

"""More convinient way of denoting the functions"""
pop_back!{T}(q::DeQue{T}) = pop!(q)

"""More convinient way of denoting the functions"""
push_front!{T}(q::DeQue{T}, x::T) = unshift!(q, x)

"""More convinient way of denoting the functions"""
push_back!{T}(q::DeQue{T}, x::T) = push!(q, x)

if !method_exists(+, Tuple{DeQue{Event}, DeQue{Event}})
    +{E}(q::DeQue{E}, p::DeQue{E}) = merge!(p, q)

@once function merge!{E}(q::DeQue{E}, p::DeQue{E})
    el = vcat(collect(q), collect(p))
    sort!(el, by=event_time)
    q = deque(E)
    for e in el
        push!(q, e)
    end
    return q
end

end

