using DataStructures
import Base.+

@once @inline front{E,O}(q::SortedSet{E,O}) = first(q)
@once back{E,O}(q::SortedSet{E,O})  = last(q)
@once @inline pop_front!{E,O}(q::SortedSet{E,O}) = pop!(q)
@once @inline pop_back!{E,O}(q::SortedSet{E,O})  = pop!(q, last(q))
@once @inline push_front!{E,O}(q::SortedSet{E,O}, e::E) = insert!(q.bt, e, nothing, true)
@once @inline push_back!{E,O}(q::SortedSet{E,O}, e::E) = insert!(q.bt, e, nothing, true)

@once @inline function merge!{E,O}(q::SortedSet{E,O}, p::SortedSet{E,O})
    for e in p
        r = (rand() - 0.5)*1e-12
        push_front!(q, Event(e.x + r, e.slope, e.offset))
    end
    q
end
