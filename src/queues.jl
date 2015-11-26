p = module_parent(current_module())
if isdefined(p, :sortedset) && typeof(p.sortedset) == Bool && p.sortedset
    info("Activating SortedSet")
    include("heap.jl")
    event_order = Base.Order.By(event_time)
    OType = typeof(event_time)
    EventQueue() = SortedSet{Event, typeof(event_order)}(event_order)
else
    info("Activating MergedArrays")
    include("depq.jl")
    EventQueue() = DePQ{Event}([], event_time)
end
