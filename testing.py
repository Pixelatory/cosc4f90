wow = [1,2,3]
tRes = [[wow[x], x, 0] for x in range(len(wow))]

print(tRes)

ok = [[[0.0,0.0], [[1,0],[0,1]]]]

ok[0].append(0)

print(ok[0][1])