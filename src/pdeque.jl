"""Prioritized double-ended queue"""
import Base.push!

type PDeQue{E}
    elements::Vector{E}
    by::Function
end

pop_front!{E}(q::PDeQue{E}) = shift!(q.elements)

pop_back!{E}(q::PDeQue{E})  = pop!(q.elements)

function push!{E}(q::PDeQue{E}, e::E)
    push!(q.elements, e)
    sort!(q.elements, by=q.by)
end
