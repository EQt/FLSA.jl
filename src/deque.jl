using DataStructures.Deque
## already done in FLSA.jl
# import DataStructures: front, back 

typealias DeQue Deque

"""More convinient way of denoting the functions"""
pop_front!{T}(q::DeQue{T}) = shift!(q)

"""More convinient way of denoting the functions"""
pop_back!{T}(q::DeQue{T}) = pop!(q)

"""More convinient way of denoting the functions"""
push_front!{T}(q::DeQue{T}, x::T) = unshift!(q, x)

"""More convinient way of denoting the functions"""
push_back!{T}(q::DeQue{T}, x::T) = push!(q, x)

+{E}(q::DeQue{E}, p::DeQue{E}) = merge!(p, q)

function merge!{E}(q::DeQue{E}, p::DeQue{E})
    el = vcat(collect(q), collect(p))
    sort!(el, by=event_time)
    q = deque(E)
    for e in el
        push!(q, e)
    end
    return q
end
