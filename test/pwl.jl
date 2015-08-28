const1 = PWL([Knot(0,1)])
@test call(const1, -1) == 1
@test call(const1, 13) == 1

constid = PWL([Knot(0,0), Knot(1,1)])
@test call(constid, -1) == -1

@test call(constid, 13) == 13

p3 = PWL([Knot(-1,-1), Knot(0,0), Knot(1,0)])
@test call(p3, -5) == -5
@test call(p3, -0.5) == -0.5
@test call(p3, 2) == 0
@test call(p3, 20) == 0
