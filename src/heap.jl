using DataStructures

@inline front{E,O}(q::SortedSet{E,O}) = first(q)
@inline back{E,O}(q::SortedSet{E,O})  = last(q)
@inline pop_front!{E,O}(q::SortedSet{E,O}) = pop!(q)
@inline pop_back!{E,O}(q::SortedSet{E,O})  = pop!(q, last(q))
@inline push_front!{E,O}(q::SortedSet{E,O}, e::E) = push!(q, e)
@inline push_back!{E,O}(q::SortedSet{E,O}, e::E) = push!(q, e)
@inline +{E,O}(q::SortedSet{E,O}, p::SortedSet{E,O}) = union!(p, q)
