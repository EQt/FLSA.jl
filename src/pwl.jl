INF = 100
EPS = 1e-9

"""
Point in a PWL
"""
type Knot
    x::Float64
    y::Float64
end

"""
Piecewise linear function.
The function is stored as a list of knots.
TODO: consider storing it as a list of (x_i, slope_i).
"""
type PWL
    knots::Vector{Knot}
end

function call(f::PWL, x::Number)
    # special cases
    if length(f.knots) == 0
        return 0
    end
    if length(f.knots) == 1
        return f.knots[1].y
    end
    i1 = findfirst(k -> k.x >= x, f.knots)
    if (i1 <= 0) # not found ==> take last ones
        p0 = f.knots[end-1]
        p1 = f.knots[end]
    elseif (i1 == 1) # first one
        p0 = f.knots[i1]
        p1 = f.knots[i1+1]
    else
        p0 = f.knots[i1-1]
        p1 = f.knots[i1]
    end
    lam = (x - p0.x) / (p1.x - p0.x)
    return (1-lam)*p0.y + lam*p1.y
end

## class PWL:
##     def __init__(self, knots):
##         """knots is list of (x_i, y_i)"""
##         self.knots = sorted(knots, key=lambda p: p[0])

##     def __call__(self, x):
##         """evaluate this function at x"""
##         try:
##             # special cases
##             if len(self.knots) == 0: return 0
##             if len(self.knots) == 1: return self.knots[0][1]
##             # find position in knots list
##             # see http://stackoverflow.com/a/2364277
##             i0 = next((i for i,p in enumerate(self.knots) if p[0] >= x)) - 1
##             i1 = i0+1
##             if i1 >= len(self.knots):
##                 i1 = len(self.knots)-1
##                 i0 = i1-1
##             p0 = self.knots[i0]
##             p1 = self.knots[i1]
##             lam = (x - p0[0]) / float(p1[0] - p0[0])
##             return (1-lam)*p0[1] + lam*p1[1]
##         except StopIteration:
##             print("x=%d, knots=%s" % (x, str(self.knots)))
##             raise
##         except IndexError:
##             print("i0=%d, i1=%d, len(knots)=%d" % (i0, i1, len(self.knots)))
##             raise
##         except ZeroDivisionError:
##             print(self)
##             print("x=%f" % x)
##             print("p0=", p0)
##             print("p1=", p1)
##             raise

function +(f::PWL, g::PWL)
    knots = vcat(f.knots, g.knots)
    sort!(knots, by=k -> k.x)
    x_old = -Inf
    off = 0
    for i = 1:length(knots)
        if (knots[i+off].x <= x_old + EPS)
            deleteat!(knots, i+off)
            off -= 1
        else
            x_old = knots[i+off].x
        end
    end
    for k in knots
        k.y = call(f, k.x) + call(g, k.x)
    end
    return PWL(knots)
end

-(f::PWL) = PWL([Knot(k.x, -k.y) for k in f.knots])
-(f::PWL, g::PWL) = f + (-g)

"Return copy(f) PWL clipped to [x_lower, x_upper]."
function clip(f::PWL, x_lower::Number, x_upper::Number)
    knots = self.knots_in_range(x_lower, x_upper)
    PWL([Knot(-INF, call(f,x_lower)), knots, Knot(+INF, call(f,x_upper))])
end


##     def find_x(self, y):
##         """find x s.t. PWL(x)=y"""
##         i0 = [i for i,p in enumerate(self.knots) if p[1] <= y][-1]
##         i1 = i0+1
##         p0, p1 = self.knots[i0], self.knots[i1]
##         lam = (y - p0[1]) / float(p1[1] - p0[1])
##         return (1-lam)*p0[0] + lam*p1[0]
