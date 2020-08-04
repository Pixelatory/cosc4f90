import random
import math
import logging
import os
import datetime
import copy
import concurrent.futures

from typing import List

from shared import aggregatedFunction, posToStrings

'''
    AMPSO for the MSA Problem
    Nick Aksamit 2020

    Acknowledgement goes towards: 
'''


def MSAAMPSO(seq, genInterval, coefLimit, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, f, w1, w2, log):
    """The AMPSO algorithm fitted for the MSA problem.

    :type seq: list of str
    :param seq: sequences to be aligned
    :type genInterval: list of float
    :param genInterval: the interval that is used within the angular modulation generation function
    :type coefLimit: list of float
    :param coefLimit: the max and min limit that all coefficients will be clamped to
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
    :type vmaxiterlimit: int
    :param vmaxiterlimit: Maximum iteration limit of clamping velocity values
    :type term: float
    :param term: termination criteria (set to float('inf') for no fitness termination)
    :type maxIter: int
    :param maxIter: maximum iteration limit (> 0)
    :type f: (list of (list of int), list of str, float, float) -> float
    :param f: fitness function (position vector, sequences, weight coefficient 1, weight coefficient 2)
    :type w1: float
    :param w1: weight coefficient for number of aligned characters
    :type w2: float
    :param w2: weight coefficient for number of leading indels used
    :type log: bool
    :param log: logging results of the MSABPSO

    :rtype: list of (list of int)
    :returns: global best position

    Initialization Process:
        Particle positions: each xi = U(coefLimit[0], coefLimit[1]), U meaning uniformly distributed random value

        Particle personal best: Same as initialized position

        Particle velocities: each vi = 0

    AMPSO:
        Topology: Star (gBest)

        PSO Type: Synchronous

        Velocity Update: v(t+1) = w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi),
        where r1 and r2 are uniformly random values between [0,1]

        Velocity clamping done until vmaxiterlimit is reached

        Position Update: x(t+1) = x(t) + v(t+1)

        Position clamping will be done for each coefficient within [coefLimit[0], coefLimit[1]]

        Angular modulation function will be within interval [genInterval[0],genInterval[1]]
    """

    # Checking for trivial errors first
    if n < 1:
        raise Exception("Swarm size cannot be < 1")
    elif len(seq) < 2:
        raise Exception("Number of sequences cannot be < 2")
    elif maxIter < 1:
        raise Exception("maxIter cannot be < 1")
    elif maxIter == float('inf') and term == float('inf'):
        raise Exception("Maximum iterations and termination fitness are both infinite!")
    elif type(coefLimit) is not list:
        raise Exception("coefLimit must be a list")
    elif len(coefLimit) < 2:
        raise Exception("coefLimit list must be at least of size 2")
    elif type(genInterval) is not list:
        raise Exception("genInterval must be a list")
    elif len(genInterval) < 2:
        raise Exception("genInterval list must be at least of size 2")

    # sort coefLimit and genInterval, just makes life easier
    if coefLimit[0] > coefLimit[1]:
        tmp = coefLimit[0]
        coefLimit[0] = coefLimit[1]
        coefLimit[1] = tmp

    if genInterval[0] > genInterval[1]:
        tmp = genInterval[0]
        genInterval[0] = genInterval[1]
        genInterval[1] = tmp

    if log:
        mkdir("MSAAMPSO logs")
        logging.basicConfig(filename=("MSAAMPSO logs/msaampso" + str(datetime.datetime.now().time()) + ".log"),
                            format='%(message)s',
                            level=logging.INFO)
        logging.info("genInterval: " + str(genInterval))
        logging.info("coefLimit: " + str(coefLimit))
        logging.info("n: " + str(n))
        logging.info("w: " + str(w))
        logging.info("c1: " + str(c1))
        logging.info("c2: " + str(c2))
        logging.info("vmax: " + str(vmax))
        logging.info("vmaxiterlimit: " + str(vmaxiterlimit))
        logging.info("term: " + str(term))
        logging.info("maxIter: " + str(maxIter))
        logging.info("f: " + str(f))
        logging.info("w1: " + str(w1))
        logging.info("w2: " + str(w2) + "\n")

    # Initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []

    def fitness(pos):
        """
        To test fitness in the AMPSO, first you use the position vector as the coefficients
        of the angular modulation formula. Then, sample random values within genInterval with
        the coefficients and use these values with the gen function. If the gen function
        returns a value > 0, the bit is 1, otherwise 0.

        :type pos: list of float
        :param pos: Position vector
        :rtype: float
        :returns: Fitness value of position vector after creating bit string with angular modulation
        """

        bitMatrix = []
        for i in range(len(seq)):
            bitMatrix.append([])
            for j in range(colLength):
                val = gen(random.uniform(genInterval[0], genInterval[1]), pos[0], pos[1], pos[2], pos[3])
                if val > 0:
                    bitMatrix[i].append(1)
                else:
                    bitMatrix[i].append(0)

        return f(bitMatrix, seq, w1, w2)

    def gen(x, a, b, c, d):
        """
        Angular Modulation Generation Function

        :type x: float
        :param x: Randomly sampled value within a range
        :type a: float
        :param a: horizontal shift coefficient
        :type b: float
        :param b: frequency coefficient
        :type c: float
        :param c: frequency coefficient
        :type d: float
        :param d: vertical shift coefficient
        :rtype: float
        """
        return math.sin(2 * math.pi * (x - a) * b * math.cos(2 * math.pi * c * (x - a))) + d

    # Just some helper variables to make code more readable
    numOfSeq = len(seq)  # number of sequences
    lSeq = {  # longest sequence (index value and length)
        "idx": 0,
        "len": len(seq[0])
    }

    # First we need to have the index and length of the longest sequence
    # logging: each sequence and its length
    for i in range(numOfSeq):
        s = len(seq[i])

        if log:
            logging.info(str(seq[i]) + " " + str(s))

        if lSeq["len"] < s:
            lSeq["idx"] = i
            lSeq["len"] = s

    if log:
        logging.info("Number of Sequences: " + str(numOfSeq) + "\n")

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength: int = math.ceil(lSeq["len"] * 1.2)

    gBestPos: List[float] = [0, 0, 0, 0]

    # Initializing the particles of swarm
    for i in range(n):
        position: List[float] = []
        velocity: List[float] = [0] * 4

        # The angular modulation coefficients are the position
        for j in range(4):
            position.append(random.uniform(coefLimit[0], coefLimit[1]))

        position.append(random.uniform(-1, 1)) # coefficient a
        position.append(random.uniform(0, 1)) # coefficient b
        position.append(random.uniform(0, 1)) # coefficient c
        position.append(random.uniform(-0.7, -0.9)) # coefficient d

        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)

        if log:
            logging.info("Particle " + str(i) + ":")
            logging.info("\tPosition: " + str(position))
            logging.info("\tPersonal Best: " + str(position))
            logging.info("\tVelocity: " + str(velocity))
            logging.info("\tFitness: " + str(fitness(position)))

        if fitness(position) > fitness(gBestPos):
            gBestPos = copy.deepcopy(position)

    if log:
        logging.info("\nGlobal best pos after particle initialization: " + str(gBestPos))
        logging.info("Global best fitness: " + str(fitness(gBestPos)) + "\n")

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter and fitness(gBestPos) < term:
        if log:
            logging.info("Iteration " + str(it))

        # Update each particle's velocity, position, and personal best
        for i in range(n):
            # r1 and r2 are ~ U (0,1)
            r1 = random.random()
            r2 = random.random()

            if log:
                logging.info("\tParticle " + str(i))
                logging.info("\t\tr1: " + str(r1))
                logging.info("\t\tr2: " + str(r2))

            # update velocity and positions in every dimension
            for j in range(4):
                pVelocities[i][j] = w * pVelocities[i][j] + \
                                    r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) + \
                                    r2 * c2 * (gBestPos[j] - pPositions[i][j])

                # velocity clamping
                if vmaxiterlimit < it:
                    if pVelocities[i][j] > vmax:
                        pVelocities[i][j] = vmax
                    elif pVelocities[i][j] < -vmax:
                        pVelocities[i][j] = -vmax

                pPositions[i][j] = pPositions[i][j] + pVelocities[i][j]

                if pPositions[i][j] > coefLimit[1]:
                    pPositions[i][j] = coefLimit[1]
                elif pPositions[i][j] < coefLimit[0]:
                    pPositions[i][j] = coefLimit[0]

            # update personal best if applicable
            if fitness(pPositions[i]) > fitness(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = copy.deepcopy(pPositions[i])

            if log:
                logging.info("\t\tPosition: " + str(pPositions[i]))
                logging.info("\t\tPersonal Best: " + str(pPersonalBests[i]))
                logging.info("\t\tVelocity: " + str(pVelocities[i]))
                logging.info("\t\tFitness: " + str(fitness(pPositions[i])))

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(n):
            if fitness(pPositions[i]) > fitness(gBestPos):  # update global best if applicable
                gBestPos = copy.deepcopy(pPositions[i])

        if log:
            logging.info("\n\tGlobal best pos: " + str(gBestPos))
            logging.info("\tGlobal best fitness: " + str(fitness(gBestPos)))

        it = it + 1

    if log:
        logging.info("\n\tFinal global best pos: " + str(gBestPos))
        logging.info("\tFinal global best fitness: " + str(fitness(gBestPos)))

    return gBestPos


def mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


# ----TESTING AREA----#
test1 = ["AGQYHECK", "AFGPWERKYV", "ASWIELKV"]
test2 = ["GAAAGTG", "CGACACTAGA", "CGCAGT"]
test3 = ["TCATGT", "GCGAT", "CGTTGT", "TCGATT", "AGCACTAG", "GAGTAGAC"]
test4 = ["DMHCMHDHMMDDMPM", "MMDCCDCCPCPCHPDPC"]
test5 = ["SCWIISRSWIWCICCRI", "WCSIWSWIWWISRICWI", "WSWWIWRCCISWCISI", "RRCCWSIRRCSRWS", "SWCRWSWSWIIRISWI"]
test6 = ["ATAHVVTAFIIWGSSGWWQFGIGVI", "IVISFVWQTIIIAGIIQFSHGAST"]
test7 = ["CDGAGIATDAWNFWAVDECVIYQIYI", "AEYGKYITDWCQLNWNCWKFTIDQGL", "GLFKLNYGDWYDVICINIQW",
         "FNADCDVYGENKETGLCAEFAENQWC", "IGGQQNLTFDLLCTIECWQYGI", "LEKQNCQNKNTTKFIIFLDDLV",
         "QIQGLYFLANGKAVVCKNKYTTN", "QFGAGFDKAEIENCQDTYCLFQGWEQK", "GFDWETLWWLIKFYEFTGTICCWNN",
         "GEDYWAGGVKIVGGICADKAEWKA"]


def testBPSOFuncWeight(seq, w1, w2):
    """Just a testing function for the BPSO on an MSA problem  (dynamic w1 and w2 values)\n
    n = 30\n
    w = 0.9\n
    c1 = c2 = 2\n
    vmax = 2\n
    vmaxiterlimit = 500\n
    term = float('inf')\n
    iter = 5000\n
    f = aggregatedFunction\n
    log = False\n

    :type seq: list of str
    :type w1: float
    :type w2: float
    :rtype: None
    """
    bestPos = []
    bestScore = 0
    sumScore = 0

    print("w1:", w1, "w2:", w2)

    e = []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in range(30):
            print("Created " + str(i))
            e.append(
                executor.submit(MSAAMPSO, seq, [-2.0, 2.0], [4.0, 4.0], 30, 0.729844, 1.49618, 1.49618, 4, 500,
                                float('inf'), 5000, aggregatedFunction, w1, w2, False))

        for future in concurrent.futures.as_completed(e):
            result = future.result()
            print("A result: " + str(result))
            score = aggregatedFunction(result, seq, w1, w2)
            sumScore += score
            print("\tScore: " + str(score))
            print("\tSum Score: " + str(sumScore))
            if score > bestScore:
                bestPos = result
                bestScore = score
            print("\tBest Pos: " + str(bestPos))
            print("\tBest Score: " + str(bestScore))

    print("Best Score:", bestScore)
    print("Average Score:", sumScore / 30)

    print("Final Result: ")
    for string in posToStrings(bestPos, seq):
        print(string)


print("test 1")
testBPSOFuncWeight(test1, 0.6, 0.4)
testBPSOFuncWeight(test1, 0.5, 0.5)
testBPSOFuncWeight(test1, 0.3, 0.7)
