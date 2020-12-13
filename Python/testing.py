#  For Math 2P91, not actually MSA stuff

#    A  B  C  D  E  F  G  H  J  K
x = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#    0  1  2  3  4  5  6  7  8  9
count = 0


def increment(x):
    global count
    i = 9
    while i >= 0 and x[i] == 9:
        x[i] = 0
        i -= 1

    if i == 0 and x[i] == 9:
        print(x, count)
        raise Exception("ok")
    else:
        x[i] += 1

    if followsRules(x):
        count += 1

def followsRules(x):
    if (x[3] * 10) + x[9] > (x[8] * 10) + x[3]:  # If DK > JD is true, then false
        return False
    elif x[4] > x[9]:  # If E > K is true, then false
        return False
    elif (x[4] * 100) + (x[3] * 10) + x[1] > (x[3] * 100) + (x[9] * 10) + x[7]:  # If EDB > DKH is true, then false
        return False
    elif 1 > x[8]:  # If 1 > J is true, then false
        return False
    elif (x[9] * 10) + x[5] > (x[3] * 10) + x[9]:  # If KF > DK is true, then false
        return False
    elif (x[4] * 1000) + (x[7] * 100) + (x[2] * 10) + x[6] > (x[4] * 1000) + (x[4] * 1000) + (x[0] * 10) + x[
        4]:  # if EHCG > EEAE is true, then false
        return False
    elif x[4] == 0:  # if E = 0 is true, then false
        return False

    # Now, checking that no two values are the same in x
    for j in range(len(x)):
        for i in range(len(x)):
            if i != j and x[i] == x[j]:
                return False

    return True


while True:
    increment(x)
