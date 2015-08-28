const1 = PWL([Knot(0,1)])
@test call(const1, -1) == 1
@test call(const1, 13) == 1

constid = PWL([Knot(0,0), Knot(1,1)])
@test call(constid, -1) == -1
@test call(constid, 13) == 13
