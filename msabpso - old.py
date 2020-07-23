import random
import math

'''
    The BPSO algorithm.

    n -> swarm size (integer > 0)
    dim -> number of dimensions (integer > 0)
    w -> inertia coefficient (float)
    c1 -> cognitive coefficient (float)
    c2 -> social coefficient (float)
    term -> termination criteria (float)
    maxIter -> maximum iteration limit (integer > 0)
    f -> the function which tests for fitness (function)

    Initialization Process:
    Particle positions -> each xi = U(0, 1), where U is uniformly distributed random value that's either 0 or 1
    Particle personal best -> Same as initialized position
    Particle velocities -> each vi = 0

    PSO Topology: Star (gBest)
    Velocity Update: w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi), where r1 and r2 are uniformly random values between [0,1]
    Probability on Velocity = p(t) = Sigmoid(v(t+1)), where v(t+1) is the new velocity vector
    Position Update: x(t+1) = {
                                1, if U(0,1) < p(t)
                                0, otherwise
                              }
'''


def MSABPSO(baseSequences, n, w, c1, c2, term, maxIter):
    # Checking for trivial errors first
    if (n < 1):
        raise Exception("Swarm size cannot be < 1")

    # initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []
    gBest = 0  # indice of the global best particle

    def f(position):
        return aggregatedFunction(position, baseSequences)

    longestSequence = 0

    # Begin with initializing the swarm
    for i in range(n):
        position = []
        velocity = []
        for j in range(len(baseSequences)):
            if len(baseSequences[longestSequence]) < len(baseSequences[j]):
                longestSequence = j

        for j in range((len(baseSequences[longestSequence]) + 1) * len(baseSequences)):
            position.append(random.randint(0, 1))
            velocity.append(0)

        position = colDashRemove(position, len(baseSequences[longestSequence]) + 1, len(baseSequences))

        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)

        if f(position) < f(pPositions[gBest]):
            gBest = i

    print("Iteration time!")
    # This is where the iterations begin
    it = 1
    while it < maxIter and f(pPositions[gBest]) > term:
        for i in range(n):  # for each particle
            r1 = random.random()
            r2 = random.random()
            for sequence in baseSequences:
                for j in range(len(sequence)):  # update the velocity and particle position using the formulas
                    pVelocities[i][j] = w * pVelocities[i][j] + \
                                        r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) + \
                                        r2 * c2 * (pPositions[gBest][j] - pPositions[i][j])
                    probability = Sigmoid(pVelocities[i][j])
                    pPositions[i][j] = 1 if random.uniform(0, 1) < probability else 0

            pPositions[i] = colDashRemove(pPositions[i], len(baseSequences[longestSequence]) + 1)

            if f(pPositions[i]) < f(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = pPositions[i]
            if f(pPositions[i]) < f(pPositions[gBest]):  # update global best if applicable
                gBest = i

        it = it + 1
    print(gBest)
    return pPositions[gBest]


'''
    colDashRemove
    
    Removes a column that consists of only dashes.
    
    x -> the particle position vector
    j -> the amount of dashes in the longest sequence within the particle position
    z -> the amount of sequences
'''
def colDashRemove(x, j, z):
    for i in range(j):
        dash = False
        if x[i] == 1:
            dash = True
            for y in range(z):
                if x[i + j * y] != 1:
                    dash = False
                    break

        if dash:
            x[i] = 0
            for y in range(z):
                x[i + j * y] = 0
    return x


def Sigmoid(x):
    return 1 / (1 + math.exp(-x))


'''
    Used for converting a position to a list of strings.
    
    position -> a vector of a particle position 
    baseSequences -> a list of sequences (strings)
    
    Long Desc:
    A position vector is one where it has a length of m * num of sequences,
    where m is the length of the longest sequence. Therefore each particle
    can have dashes
'''


def posToStrings(position, baseSequences):
    result = []
    for i in range(len(baseSequences)):
        start = i * (len(position) / len(baseSequences))
        end = start + len(baseSequences[i]) + 1
        bits = position[int(start):int(end)]
        string = ""
        offset = 0
        for bit in bits:
            if bit == 1:
                string = string + "-"

            if offset < len(baseSequences[i]):
                string = string + baseSequences[i][offset]
                offset = offset + 1
        result.append(string)
    return result


def numOfAlignedChars(strings):
    if len(strings) < 1:
        raise Exception("There's no strings in the num of aligned chars function")
    elif len(strings) == 1:
        print("Warning: only 1 string in the numOfAlignedChars function")
        return 0

    result = 0
    longestString = 0
    for i in range(len(strings)):
        if len(strings[longestString]) < len(strings[i]):
            longestString = i

    for i in range(len(strings[longestString])):
        for j in range(len(strings)):
            for string in strings[j + 1:]:
                try:
                    if strings[j][i] == string[i] and strings[j][i] != "-":
                        result = result + 1
                        # print(i, strings[j], string, "true")
                except IndexError:
                    # print(i, strings[j], string, "false")
                    pass
    return result


def numOfSpaces(strings):
    result = 0
    for string in strings:
        for char in string:
            if char == "-":
                result = result + 1
    return result


def aggregatedFunction(position, baseSequences):
    strings = posToStrings(position, baseSequences)
    return numOfAlignedChars(strings) - numOfSpaces(strings)

'''
baseSequences = ["AABGT", "AABGT", "ALWI"]

bestPos = MSABPSO(baseSequences, 15, 0.7, 1.4, 1.4, 500, 5000)
bestStrings = posToStrings(bestPos, baseSequences)

print("Best score: " + str(aggregatedFunction(bestPos, baseSequences)))
for string in bestStrings:
    print(string)
'''
