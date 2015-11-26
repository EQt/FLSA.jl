using DataStructures
import Base.merge!

@inline front{E,O}(q::SortedMultiDict{E,Void,O}) = first(q)
@inline back{E,O}(q::SortedMultiDict{E,Void,O})  = last(q)
@inline pop_front!{E,O}(q::SortedMultiDict{E,Void,O}) = pop!(q)
@inline pop_back!{E,O}(q::SortedMultiDict{E,Void,O})  = pop!(q, last(q))
@inline push_front!{E,O}(q::SortedMultiDict{E,Void,O}, e::E) = insert!(q, e, nothing)
@inline push_back!{E,O}(q::SortedMultiDict{E,Void,O}, e::E) = insert!(q, e, nothing)
@inline function +{E,O}(p::SortedMultiDict{E,Void,O}, q::SortedMultiDict{E,Void,O})
    merge!(p, q)
    p
end
