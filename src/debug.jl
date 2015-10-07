"""Debugging functionality of FLSA"""

p = module_parent(current_module())
if isdefined(p, :debug) && typeof(p.debug) == Bool && p.debug == true
    macro debug(msg)
        :(println("DEBUG: " * $msg))
    end
elseif isdefined(p, :debug) && typeof(p.debug) == Function
    info("Activating DEBUGGING in FLSA")
    macro debug(msg)
        :(p.debug($msg))
    end
else
    # info("Deactivating debugging in FLSA")
    macro debug(msg)
        # do nothing
    end
end
