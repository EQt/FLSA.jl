"""Debugging functionality of FLSA"""

p = module_parent(current_module())
if true || isdefined(p, :debug) && typeof(p.debug) == Function
    info("Activating DEBUGGING in FLSA")
    macro debug(msg)
        :(println("DEBUG: " * $msg))
    end
else
    info("Deactivating debugging in FLSA")
    macro debug(msg)
        # do nothing
    end
end
