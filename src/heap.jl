using DataStructures

@inline front{E,O}(q::SortedSet{E,O}) = first(q)
@inline back{E,O}(q::SortedSet{E,O})  = last(q)
@inline pop_front!{E,O}(q::SortedSet{E,O}) = pop!(q)
@inline pop_back!{E,O}(q::SortedSet{E,O})  = pop!(q, last(q))
@inline push_front!{E,O}(q::SortedSet{E,O}, e::E) = insert!(q.bt, e, nothing, true)
@inline push_back!{E,O}(q::SortedSet{E,O}, e::E) = insert!(q.bt, e, nothing, true)
@inline function +{E,O}(q::SortedSet{E,O}, p::SortedSet{E,O})
    for e in p
        push_front!(q, e)
    end
    q
end

