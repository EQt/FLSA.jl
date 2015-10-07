"""More convinient way of denoting the functions"""
using DataStructures.Deque

pop_front!{T}(q::Deque{T}) = shift!(q)
pop_back!{T}(q::Deque{T}) = pop!(q)
push_front!{T}(q::Deque{T}, x::T) = unshift!(q, x)
push_back!{T}(q::Deque{T}, x::T) = push!(q, x)
