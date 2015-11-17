import Base:+, -

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
    PWL(x::Float64, y::Float64; slope::Float64 = 1.0) = new([Knot(x-INF,y-INF), Knot(x+INF, y+INF)])
    PWL(ks::Vector{Knot}) = new(ks)
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

function +(f::PWL, g::PWL)
    knots = deepcopy([f.knots; g.knots])
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

"""Return clipped f to x-range [x_lower, x_upper] (unlike clip in text)"""
function clip_x(f::PWL, x_lower::Number, x_upper::Number)
    i0 = findfirst(k -> k.x >= x_lower, f.knots)
    i1 = findlast( k -> k.x <= x_upper, f.knots)
    knots = [Knot(x_lower, call(f, x_lower));
             f.knots[i0:i1];
             Knot(x_upper, call(f, x_upper))]
    PWL([Knot(-INF, call(f,x_lower)); knots; Knot(+INF, call(f,x_upper))])
end


"""Return clipped f to range [lower, upper]"""
clip(f::PWL, lower::Number, upper::Number) = clip(f, find_x(f, lower), find_x(f, upper))


"""find x s.t. PWL(x)=y (assumes that f is strictly increasing)"""
function find_x(f::PWL, y::Number)
    i0 = findlast(k -> k.y <= y, f.knots)
    if (i0 == findlast(f.knots)) i0 -= 1 end
    i1 = i0+1
    p0, p1 = f.knots[i0], f.knots[i1]
    lam = (y - p0.y) / (p1.y - p0.y)
    return (1-lam)*p0.x + lam*p1.x
end
