"""Prioritized double-ended queue"""
import Base.push!
import Base.start
import Base.next
import Base.done
import Base.length
import Base.isempty

type PDeQue{E}
    elements::Vector{E}
    by::Function
end

front{E}(q::PDeQue{E}) = q.elements[1]

back{E}(q::PDeQue{E})  = q.elements[end]

pop_front!{E}(q::PDeQue{E}) = shift!(q.elements)

pop_back!{E}(q::PDeQue{E})  = pop!(q.elements)

function push!{E}(q::PDeQue{E}, e::E)
    push!(q.elements, e)
    sort!(q.elements, by=q.by)
end

length{E}(q::PDeQue) = length(q.elements)
isempty{E}(q::PDeQue) = length(q) == 0

# for convinience (map, enumerate, …)
start{E}(q::PDeQue{E}) = start(q.elements)
next{E}(q::PDeQue{E}, s) = next(q.elements, s)
done{E}(q::PDeQue{E}, s) = done(q.elements, s)
