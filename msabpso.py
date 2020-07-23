import random
import math


def MSABPSO(seq, n, w, c1, c2, vmax, term, maxIter, f, w1, w2):
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
    :type vmax: float
    :param vmax: maximum velocity value (clamping) (set to float('inf') for no clamping)
    :type term: float
    :param term: termination criteria (set to float('inf') for no fitness termination)
    :type maxIter: int
    :param maxIter: maximum iteration limit (> 0)
    :type f: ((list of (list of int), list of str, float, float) -> float)
    :param f: fitness function
    :type w1: float
    :param w1: weight coefficient for number of aligned characters
    :type w2: float
    :param w2: weight coefficient for number of leading indels used

    :rtype: list of (list of int)
    :returns: global best position

    Initialization Process:
        Particle positions: each xi = U(0, 1), where U is uniformly distributed random value that's either 0 or 1

        Particle personal best: Same as initialized position

        Particle velocities: each vi = 0

    BPSO:
        Topology: Star (gBest)

        PSO Type: Synchronous

        Velocity Update: v(t+1) = w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi), where r1 and r2 are uniformly random values between [0,1]

        Probability on Velocity: p(t+1) = Sigmoid(v(t+1)), where v(t+1) is the new velocity vector

        Position Update: x(t+1) = { 1, if U(0,1) < p(t+1); 0, otherwise }
    """

    # Checking for trivial errors first
    if n < 1:
        raise Exception("Swarm size cannot be < 1")
    elif len(seq) < 2:
        raise Exception("Number of sequences < 2")
    elif maxIter < 1:
        raise Exception("maxIter cannot be < 1")
    elif maxIter == float('inf') and term == float('inf'):
        raise Exception("Maximum iterations and termination fitness are both infinite!")

    # Initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []
    gBest = 0  # indice of the global best particle

    def fitness(pos):
        """
        :type pos: list of (list of int)
        :param pos: Position vector
        :rtype: float
        :returns: Fitness value of position vector
        """
        return f(pos, seq, w1, w2)

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
        for j in range(numOfSeq):
            position.append([0] * colLength)
            velocity.append([0] * colLength)
            for x in range(colLength - len(seq[j])):
                while True:
                    randNum = random.randint(0, colLength - 1)
                    if position[j][randNum] != 1:
                        position[j][randNum] = 1
                        break

        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)

        if fitness(position) > fitness(pPositions[gBest]):
            gBest = i

    # This is where the iterations begin

    it = 0  # iteration count
    while it < maxIter and fitness(pPositions[gBest]) > term:

        # Update each particle's velocity, position, and personal best
        for i in range(n):
            # r1 and r2 are ~ U (0,1)
            r1 = random.random()
            r2 = random.random()

            # update velocity and positions in every dimension
            for j in range(numOfSeq):
                for x in range(colLength):
                    pVelocities[i][j][x] = w * pVelocities[i][j][x] + \
                                           r1 * c1 * (pPersonalBests[i][j][x] - pPositions[i][j][x]) + \
                                           r2 * c2 * (pPositions[gBest][j][x] - pPositions[i][j][x])

                    # velocity clamping
                    if pVelocities[i][j][x] > vmax:
                        pVelocities[i][j][x] = vmax
                    elif pVelocities[i][j][x] < -vmax:
                        pVelocities[i][j][x] = -vmax

                    probability = Sigmoid(pVelocities[i][j][x])
                    pPositions[i][j][x] = 1 if random.uniform(0, 1) < probability else 0

            # update personal best if applicable
            if fitness(pPositions[i]) > fitness(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = pPositions[i]

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(n):
            if fitness(pPositions[i]) > fitness(pPositions[gBest]):  # update global best if applicable
                gBest = i

        it = it + 1
    return pPositions[gBest]

def Sigmoid(x):
    """The classic sigmoid function.

    :type x: int
    :rtype: float
    """
    return 1 / (1 + math.exp(-x))

def posToStrings(position, seq):
    """Converts a list of sequences into a list of strings with indels, according to the position vector given.

    :type position: list of (list of int)
    :type seq: list of str
    :rtype: list of str
    """
    result = []
    i = 0
    for bitlist in position:
        j = 0
        result.append("")
        for bit in bitlist:
            if bit == 0 and j < len(seq[i]):
                result[len(result) - 1] = result[len(result) - 1] + seq[i][j]
                j = j + 1
            else:
                result[len(result) - 1] = result[len(result) - 1] + "-"
        i = i + 1
    return result


def numOfAlignedChars(strings):
    """Counts the number of aligned characters in a list of strings.

    :type strings: list of str
    :rtype: int
    """
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

def aggregatedFunction(position, seq, w1, w2):
    """A maximization aggregated fitness function that follows the following formula:

    f(x) = w1 * numOfAlignedChars(x) + w2 * (nMax - nI),\n
    where nMax is the number of total indels,\n
    and nI is the number of indels in-between characters.

    Note: if the position vector is invalid, then -inf is returned


    :param position: position vector
    :type position: list of (list of int)
    :param seq: sequences to be aligned
    :type seq: list of str
    :param w1: weight coefficient for number of aligned characters
    :type w1: float
    :param w2: weight coefficient for number of leading indels used
    :type w2: float
    :rtype: float
    :return: fitness value
    """

    strings = posToStrings(position, seq)

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
    for i in range(len(seq)):
        tmp = 0
        hitFirstChar = False
        for bit in position[i]:
            if bit == 0:
                tmp = tmp + 1
                hitFirstChar = True
            elif bit == 1 and hitFirstChar and tmp < len(seq[i]):
                nI = nI + 1  # an indel was found in-between characters

        if tmp != len(seq[i]):
            return float('-inf')  # return a very small number, this position is invalid

    return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI))

# unused as of right now
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



#----TESTING AREA----#
baseSequences = ["CDGAGIATDAWNFWAVDECVIYQIYI", "AEYGKYITDWCQLNWNCWKFTIDQGL", "GLFKLNYGDWYDVICINIQW", "FNADCDVYGENKETGLCAEFAENQWC", "IGGQQNLTFDLLCTIECWQYGI", "LEKQNCQNKNTTKFIIFLDDLV", "QIQGLYFLANGKAVVCKNKYTTN", "QFGAGFDKAEIENCQDTYCLFQGWEQK", "GFDWETLWWLIKFYEFTGTICCWNN", "GEDYWAGGVKIVGGICADKAEWKA"]

def testBPSOFuncWeight(seq, f, w1, w2):
    bestPos = []
    bestScore = 0

    print("w1:", w1, "w2:", w2)

    for i in range(30):
        pos = MSABPSO(seq, 30, 0.9, 2, 2, 2, float('inf'), 5000, f, 0.5, 0.5)
        score = f(pos,seq,w1,w2)
        if score > bestScore:
            bestScore = score
            bestPos = pos

    print("Best Score:", bestScore)
    for string in posToStrings(bestPos,seq):
        print(string)

testBPSOFuncWeight(baseSequences, aggregatedFunction, 0.6, 0.4)
testBPSOFuncWeight(baseSequences, aggregatedFunction, 0.5, 0.5)
testBPSOFuncWeight(baseSequences, aggregatedFunction, 0.3, 0.7)