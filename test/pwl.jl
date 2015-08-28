const1 = PWL(Knot(0,1))
@test call(const1, -1) == 1
@test call(const1, 13) == 1
