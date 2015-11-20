using DataStructures.Deque
import DataStructures: front, back

typealias DeQue Deque

"""More convinient way of denoting the functions"""
pop_front!{T}(q::DeQue{T}) = shift!(q)

"""More convinient way of denoting the functions"""
pop_back!{T}(q::DeQue{T}) = pop!(q)

"""More convinient way of denoting the functions"""
push_front!{T}(q::DeQue{T}, x::T) = unshift!(q, x)

"""More convinient way of denoting the functions"""
push_back!{T}(q::DeQue{T}, x::T) = push!(q, x)
