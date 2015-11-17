using DataStructures.Deque

typealias DeQue Deque

"""More convinient way of denoting the functions"""
pop_front!{T}(q::DeQue{T}) = shift!(q)

pop_back!{T}(q::DeQue{T}) = pop!(q)

push_front!{T}(q::DeQue{T}, x::T) = unshift!(q, x)

push_back!{T}(q::DeQue{T}, x::T) = push!(q, x)
