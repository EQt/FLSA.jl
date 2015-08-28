INF = 100
EPS = 1e-9

type Knot
    x::Float64
    y::Float64
end

type PWL
    """
    Piecewise linear function.
    The function is stored as a list of knots.
    TODO: consider storing it as a list of (x_i, slope_i).
    """
    knots::Vector{Knot}
end

using Base
function call(f::PWL, x::Number)
    # special cases
    if length(f.knots) == 0
        return 0
    end
    if length(f.knots) == 1
        return f.knots[1].y
    end
end

testf = PWL([Knot(1,1), Knot(2,1)])

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

##     def __add__(self, other):
##         """add to PWLs"""
##         if isinstance(other, Number):
##             return PWL([(x, y+other) for x, y in self.knots])
##         if not isinstance(other, PWL):
##             raise TypeError('Unsupported type ' + str(type(other)))
##         xs = sorted(map(lambda p: p[0], self.knots + other.knots))
##         xs_unique, x_old = [], float('-inf')
##         for x in xs:
##             if x > x_old + EPS:
##                 xs_unique.append(x)
##                 x_old = x
##         rsum = lambda x: self(x) + other(x)
##         knots = [(x, rsum(x)) for x in xs_unique]
##         return PWL(knots)

##     def __radd__(self, other):
##         if isinstance(other, Number):
##             return PWL([(x, y+other) for x, y in self.knots])
##         else:
##             raise RuntimeError('this should not happen')

##     def __neg__(self):
##         """mirror at x axis"""
##         return PWL([(x, -y) for x,y in self.knots])

##     def __sub__(self, other):
##         """substract to PWLs"""
##         return self + (-other)

##     def __str__(self):
##         return "\n".join(map(lambda p: "(% 8.2f, % 8.2f)" % p, self.knots))

##     def clip(self, x_lower, x_upper):
##         """Return copy(self) PWL clipped to [x_lower, x_upper]."""
##         knots = self.knots_in_range(x_lower, x_upper)
##         return PWL([(-INF, self(x_lower))] + knots + [(+INF, self(x_upper))])


##     def plot(self, xmin=-INF, xmax=+INF):
##         """plot using matplotlib"""
##         import matplotlib.pyplot as plt
##         knots = [k for k in self.knots if xmin <= k[0] and k[0] <= xmax]
##         knots = [(xmin, self(xmin))] + knots + [(xmax, self(xmax))]
##         plt.plot([k[0] for k in knots], [k[1] for k in knots], '-o')


##     def find_x(self, y):
##         """find x s.t. PWL(x)=y"""
##         i0 = [i for i,p in enumerate(self.knots) if p[1] <= y][-1]
##         i1 = i0+1
##         p0, p1 = self.knots[i0], self.knots[i1]
##         lam = (y - p0[1]) / float(p1[1] - p0[1])
##         return (1-lam)*p0[0] + lam*p1[0]


##     def knots_in_range(self, xmin, xmax):
##         """return knots within the range [xmin, xmax]"""
##         mid = [k for k in self.knots if xmin <= k[0] and k[0] <= xmax]
##         return [(xmin, self(xmin))] + mid + [(xmax, self(xmax))]
