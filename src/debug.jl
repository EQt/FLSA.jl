"""Debugging functionality of FLSA"""

p = module_parent(current_module())
if isdefined(p, :debug) && typeof(p.debug) == Bool && p.debug == true
    info("Activating DEBUGGING in FLSA")
    macro debug(msg)
        :(println("DEBUG: " * $msg))
    end
elseif isdefined(p, :debug) && typeof(p.debug) == Function
    info("Activating DEBUGGING in FLSA for $(p.debug.env.defs)")
    macro debug(msg)
        :(p.debug($msg))
    end
else
    # info("Deactivating debugging in FLSA (p = $p)")
    macro debug(msg)
        # do nothing
    end
    ## import Base.@assert
    ## macro assert(ex, msg...)
    ##     # do nothing
    ## end
end


macro val(v)
    s = "$v="
     :($s * string($v))
end


macro once(expr)
    return quote
        if !isdefined(:LastMain) || !isdefined(LastMain, :FLSA)
            info("Defining first time")
            $expr
        else
            info("Already defined")
        end
    end
end

