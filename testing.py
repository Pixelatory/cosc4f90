
def algorithm(tStopper):
    t=1
    s=t
    n=1

    loops = 0

    while True:
        loops += 1

        n = n + 1
        t *= (n-1) / (2*n-1)
        s += t

        print(t)

        if t < tStopper:
            break
    print()
    print(2*s, loops)

algorithm(10**-4)
