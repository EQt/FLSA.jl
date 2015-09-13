"""Debugging functionality of FLSA"""

if isdefined(:debug)
    macro debug(msg)
        :(info($msg))
    end
else
    macro debug(msg)
    end
end
