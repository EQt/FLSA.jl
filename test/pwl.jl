const1 = PWL([Knot(0,1)])
@test call(const1, -1) == 1
@test call(const1, 13) == 1

constid = PWL([Knot(0,0); Knot(1,1)])
@test call(constid, -1) == -1

@test call(constid, 13) == 13

p3 = PWL([Knot(-1,-1); Knot(0,0); Knot(1,0)])
@test call(p3, -5) == -5
@test call(p3, -0.5) == -0.5
@test call(p3, 2) == 0
@test call(p3, 20) == 0

p4 = p3 + constid
@test call(p4, -5) == -10.0
@test call(p4, -0.5) == -1.0
@test call(p4, 2) == 2
@test call(p4, 20) == 20

# p3 should stay the same...
@test call(p3, -5) == -5
@test call(p3, -0.5) == -0.5
@test call(p3, 2) == 0
@test call(p3, 20) == 0


@test find_x(p3,   -1) ==   -1
@test find_x(p3, -0.5) == -0.5
@test find_x(p4,   20) ==   20
@test find_x(p4, -1.0) == -0.5

con = clip(constid, -1, 1)
@test call(con, -2) == -1
@test call(con, 0.5) == 0.5
@test call(con, 13) == 1
