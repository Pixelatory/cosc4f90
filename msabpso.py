import random
import math


def MSABPSO(seq, n, w, c1, c2, term, maxIter, f):
    """The BPSO algorithm fitted for the MSA problem.

    :type seq: list of str
    :param seq: sequences to be aligned
    :type n: int
    :param n: swarm size (> 0)
    :type w: float
    :param w: inertia coefficient
    :type c1: float
    :param c1: cognitive coefficient
    :type c2: float
    :param c2: social coefficient
    :type term: float
    :param term: termination criteria
    :type maxIter: int
    :param maxIter: maximum iteration limit (> 0)
    :type f: function
    :param f: fitness function

    :rtype: list[int]
    :returns: The global best

    Initialization Process:
        Particle positions -> each xi = U(0, 1), where U is uniformly distributed random value that's either 0 or 1

        Particle personal best -> Same as initialized position

        Particle velocities -> each vi = 0

    BPSO:
        Topology -> Star (gBest)

        Velocity Update -> v(t+1) = w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi), where r1 and r2 are uniformly random values between [0,1]

        Probability on Velocity -> p(t) = Sigmoid(v(t+1)), where v(t+1) is the new velocity vector

        Position Update -> x(t+1) = { 1, if U(0,1) < p(t); 0, otherwise }
    """

    # Checking for trivial errors first
    if n < 1:
        raise Exception("Swarm size cannot be < 1")
    elif len(seq) < 2:
        raise Exception("Number of sequences < 2")
    elif maxIter < 1:
        raise Exception("maxIter cannot be < 1")

    # Initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []
    gBest = 0  # indice of the global best particle

    def fitness(pos):
        """
        :type pos: list of int
        :param pos: Position vector
        :rtype: float
        :returns: Fitness value of position vector
        """
        return f(pos, seq)

    # Just some helper variables to make code more readable
    numOfSeq = len(seq)  # number of sequences
    lSeq = {  # longest sequence (index value and length)
        "idx": 0,
        "len": len(seq[0])
    }

    # First we need to have the index and length of the longest sequence
    for i in range(numOfSeq):
        s = len(seq[i])
        if lSeq["len"] < s:
            lSeq["idx"] = i
            lSeq["len"] = s

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength = math.ceil(lSeq["len"] * 1.2)

    # Initializing the particles of swarm
    for i in range(n):
        position = []
        velocity = []

        # This initiates the position and velocity, which are
        # both a matrix of size [numOfSequences x colLength].
        # Additionally it randomly puts indels in a random spot,
        # but only just the right amount of them so the initial
        # position is feasible.
        for s in seq:
            position.append([0] * colLength)
            velocity.append([0] * colLength)
            for x in range(colLength - len(s)):
                while True:
                    randNum = random.randint(0, colLength)
                    if position[randNum] != 1:
                        position[randNum] = 1
                        break

        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)

        if fitness(position) > fitness(pPositions[gBest]):
            gBest = i

    # This is where the iterations begin

    print("Swarm done initialization, now iteration time!")
    it = 0  # iteration count
    while it < maxIter and fitness(pPositions[gBest]) > term:

        for i in range(n):  # for each particle
            # generating the random stochastic scalars
            r1 = random.random()
            r2 = random.random()

            for j in range(numOfSeq):  # updating velocity and position (using probability)
                for x in range(colLength):
                    pVelocities[i][j][x] = w * pVelocities[i][j][x] + \
                                           r1 * c1 * (pPersonalBests[i][j][x] - pPositions[i][j][x]) + \
                                           r2 * c2 * (pPositions[gBest][j][x] - pPositions[i][j][x])
                    probability = Sigmoid(pVelocities[i][j][x])
                    pPositions[i][j][x] = 1 if random.uniform(0, 1) < probability else 0

            if fitness(pPositions[i]) > fitness(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = pPositions[i]
            if fitness(pPositions[i]) > fitness(pPositions[gBest]):  # update global best if applicable
                gBest = i

        it = it + 1
    return pPositions[gBest]


'''
    Sigmoid
    
    The classic sigmoid function.
'''


def Sigmoid(x):
    return 1 / (1 + math.exp(-x))


'''
    Used for converting a position to a list of strings.
    
    position -> a vector of a particle position 
    baseSequences -> a list of sequences (strings)
'''


def posToStrings(position, baseSequences):
    result = []
    i = 0
    for bitlist in position:
        j = 0
        result.append("")
        for bit in bitlist:
            if bit == 0 and j < len(baseSequences[i]):
                result[len(result) - 1] = result[len(result) - 1] + baseSequences[i][j]
                j = j + 1
            else:
                result[len(result) - 1] = result[len(result) - 1] + "-"
        i = i + 1
    return result


def numOfAlignedChars(strings):
    if len(strings) < 1:
        raise Exception("There's no strings in the num of aligned chars function")
    elif len(strings) == 1:
        print("Warning: only 1 string in the numOfAlignedChars function")
        return 0

    result = 0
    charList = []
    for i in range(len(strings[0])):
        charList.append({})
        for string in strings:
            if string[i] != "-":
                if string[i] in charList[i]:
                    charList[i][string[i]] = charList[i][string[i]] + 1
                else:
                    charList[i][string[i]] = 1

    for dict in charList:
        for v in dict.values():
            if v > 1:
                result = result + v

    return result


'''
    aggregatedFunction - to be maximized
    
    position -> position vector
    baseSequences -> list of the base sequences used for position vector
    w1 -> weight factor for number of aligned characters
    w2 -> weight factor for number of leading indels
    
    If the position given is invalid for the MSA, then -infinity will be
    returned.
'''


def aggregatedFunction(position, baseSequences, w1, w2):
    strings = posToStrings(position, baseSequences)

    nMax = 0  # total number of indels

    for bitlist in position:
        for bit in bitlist:
            if bit == 1:
                nMax = nMax + 1

    nI = 0  # number of indels inbetween chars

    # The small procedure below counts for nI,
    # and also eliminates infeasible positions.
    # If the number of 0 bits is more than the
    # amount of characters for the sequence, then
    # it's invalid. (return -infinity)
    for i in range(len(baseSequences)):
        tmp = 0
        hitFirstChar = False
        for bit in position[i]:
            if bit == 0:
                tmp = tmp + 1
                hitFirstChar = True
            elif bit == 1 and hitFirstChar and tmp < len(baseSequences[i]):
                nI = nI + 1  # an indel was found in-between characters

        if tmp != len(baseSequences[i]):
            return float('-inf')  # return a very small number, this position is invalid

    return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI))


'''
    UNUSED AS OF RN
    colDashRemove

    Removes a column that consists of only dashes.

    x -> the position vector
    y -> the sequences list
'''


def colDashRemove(x, y):
    for i in range(len(x[0])):
        dash = False
        if x[0][i] == 1:
            dash = True
            for bitlist in x:
                if bitlist[i] != 1:
                    dash = False
                    break

        if dash:
            print(str(i) + " is dash")
            for bitlist in x:
                bitlist[i] = 0
    return x


pos = [
    [0, 0, 1, 0, 1, 1, 0],
    [1, 1, 0, 0, 0, 1, 1]
]

seq = ["AFOW", "WOS"]
for string in posToStrings(pos, seq):
    print(string)

print(aggregatedFunction(pos, seq, 1, 1))

print(type(>))

baseSequences = ["AABGT", "AABGT", "ALWI"]

bestPos = MSABPSO(baseSequences, 15, 0.7, 1.4, 1.4, 500, 5000)
bestStrings = posToStrings(bestPos, baseSequences)

print("Best score: " + str(aggregatedFunction(bestPos, baseSequences)))
for string in bestStrings:
    print(string)
