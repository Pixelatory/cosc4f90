import random
import math
import datetime
import concurrent.futures
import logging

from copy import deepcopy
from statistics import stdev
from operator import lt, gt
from typing import List
from util import aggregatedFunction, getLongestSeqDict, numOfAlignedChars, bitsToStrings, numOfInsertedIndels, \
    genBitMatrix, test1, test2, test3, test4, test5, test6, test7, infeasible

"""
    AMPSO for the MSA Problem
    Nick Aksamit 2020

    Acknowledgement goes towards: 
"""


def MSAAMPSO(seq, genInterval, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, f, w1, w2, ops):
    """The AMPSO algorithm fitted for the MSA problem.

    Initialization Process:
        Particle positions: each xi = TODO replace this

        Particle personal best: Same as initialized position

        Particle velocities: each vi = 0

    AMPSO:
        Topology: Star (gBest)

        PSO Type: Synchronous

        Velocity Update: v(t+1) = w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi),
        where r1 and r2 are uniformly random values between [0,1]

        Velocity clamping done until vmaxiterlimit is reached

        Position Update: x(t+1) = x(t) + v(t+1)

        Angular modulation function will be within interval [genInterval[0], genInterval[1]]


        :type seq: List[str]
        :param seq: sequences to be aligned
        :type genInterval: List[float]
        :param genInterval: the interval that is used within the angular modulation generation function
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
        :type f: (List[List[int]], List[str], float, float, bool, List[(int, int) -> bool]) -> float
        :param f: fitness function (bitmatrix, sequences, weight coefficient 1, weight coefficient 2, checkInfeasability,
        operator)

        :type w1: float
        :param w1: weight coefficient for number of aligned characters
        :type w2: float
        :param w2: weight coefficient for number of leading indels used
        :type ops: List[(float, float) -> bool]
        :param ops: operators that will check for infeasibility (See util.py -> infeasible)
        :rtype: (List[float], List[List[int]], int)
        :returns: (global best position, global best bitmatrix, numOfInfeasibleSols)
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
    elif type(genInterval) is not list:
        raise Exception("genInterval must be a list")
    elif len(genInterval) < 2:
        raise Exception("genInterval list must be at least of size 2")

    # sort genInterval, just makes life easier
    if genInterval[0] > genInterval[1]:
        tmp = genInterval[0]
        genInterval[0] = genInterval[1]
        genInterval[1] = tmp

    # Initialize the data containers
    pPositions: List[List[float]] = []
    pPersonalBests: List[List[float]] = []
    pVelocities: List[List[float]] = []
    pBitStrings: List[List[List[int]]] = []
    numOfInfeasibleSols: int = 0

    def fitness(bitmatrix):
        """
        To test fitness in the AMPSO, first you use the position vector as the coefficients
        of the angular modulation formula. Then, sample random values within genInterval with
        the coefficients and use these values with the gen function. If the gen function
        returns a value > 0, the bit is 1, otherwise 0.

        :type bitmatrix: List[List[int]]
        :param bitmatrix: Two-dimensional binary matrix
        :rtype: float
        :returns: Fitness value of bit string
        """

        return f(bitmatrix, seq, w1, w2, True, ops)

    lSeq = getLongestSeqDict(seq)  # Longest sequence value dictionary

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength: int = math.ceil(lSeq["len"] * 1.2)

    # Ensure we start with a valid gBest particle
    while True:
        gBest = {
            "pos": [random.uniform(-0.3, 0.3), random.uniform(0.5, 10), random.uniform(0.5, 10),
                    random.uniform(-0.9, -0.5),
                    random.uniform(-0.9, -0.5)]}
        gBest["bitstring"] = genBitMatrix(gBest["pos"], seq, colLength, genInterval)

        if fitness(gBest["bitstring"]) != -float('inf'):
            break

    # Initializing the particles of swarm
    for i in range(n):
        position: List[float] = []
        velocity: List[float] = [0] * 4

        position.append(random.uniform(-0.3, 0.3))  # coefficient a
        position.append(random.uniform(0.5, 10))  # coefficient b
        position.append(random.uniform(0.5, 10))  # coefficient c
        position.append(random.uniform(-0.9, -0.5))  # coefficient d

        bitstring = genBitMatrix(position, seq, colLength, genInterval)

        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)
        pBitStrings.append(bitstring)

        tmpScore = fitness(bitstring)

        # Infeasible solution, increment numOfInfeasibleSols count and don't test for gBest
        if tmpScore == -float('inf'):
            numOfInfeasibleSols += 1
            continue

        if tmpScore > fitness(gBest["bitstring"]):
            gBest["pos"] = deepcopy(position)
            gBest["bitstring"] = deepcopy(bitstring)

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter and fitness(gBest["bitstring"]) < term:

        # Update each particle's velocity, position, and personal best
        for i in range(n):
            # r1 and r2 are ~ U (0,1)
            r1 = random.random()
            r2 = random.random()

            # update velocity and positions in every dimension
            for j in range(4):
                pVelocities[i][j] = w * pVelocities[i][j] + \
                                    r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) + \
                                    r2 * c2 * (gBest["pos"][j] - pPositions[i][j])

                # velocity clamping
                if vmaxiterlimit < it:
                    if pVelocities[i][j] > vmax:
                        pVelocities[i][j] = vmax
                    elif pVelocities[i][j] < -vmax:
                        pVelocities[i][j] = -vmax

                pPositions[i][j] += pVelocities[i][j]

            bitstring = genBitMatrix(pPositions[i], seq, colLength, genInterval)

            # update personal best if applicable
            tmpScore = fitness(bitstring)

            # Infeasible solution, increment numOfInfeasibleSols count and don't test for personal best
            if tmpScore == -float('inf'):
                numOfInfeasibleSols += 1
                continue

            # Feasible sol, so update personal best if applicable
            if tmpScore > fitness(pBitStrings[i]):
                pPersonalBests[i] = deepcopy(pPositions[i])
                pBitStrings[i] = bitstring

        # update the global best after all positions were changed (synchronous PSO)
        # Notice: no need to check for infeasibility here, cause if it was, wouldn't be added to personal best anyways
        for i in range(n):
            if fitness(pBitStrings[i]) > fitness(gBest["bitstring"]):  # update global best if applicable
                gBest["pos"] = deepcopy(pPositions[i])
                gBest["bitstring"] = deepcopy(pBitStrings[i])

        it = it + 1

    return gBest["pos"], gBest["bitstring"], numOfInfeasibleSols


# ----TESTING AREA----#


def testAMPSOFuncWeight(seq, w1, w2):
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
    sumInserted = 0
    sumAligned = 0
    bestInserted = 0
    bestAligned = 0
    bestBitString = []
    scores = []
    inserts = []
    aligns = []

    # Just logging stuff here and printing to screen
    logging.info("Started " + str(datetime.datetime.now().time()))
    logging.info("w1: " + str(w1) + " w2: " + str(w2))
    print("Started " + str(datetime.datetime.now().time()))
    print("w1:", w1, "w2:", w2)

    e = []

    # Multi-threading all 30 runs of the testing
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in range(30):
            e.append(
                executor.submit(MSAAMPSO, seq, [-100, 100], 30, 0.7, 2,
                                2, float('inf'), 500, float('inf'), 5000, aggregatedFunction, w1, w2, [lt]))

        for future in concurrent.futures.as_completed(e):
            result = future.result()
            pos = result[0]
            matrix = result[1]
            numOfInfeasibleSols = result[2]

            # Just logging and printing to screen
            logging.info("A result: " + str(pos))
            logging.info("Infeasible Sol: " + str((numOfInfeasibleSols / (30 * 5000)) * 100) + "%")
            print("A result: " + str(pos))
            print("Infeasible Sol: " + str((numOfInfeasibleSols / (30 * 5000)) * 100) + "%")

            score = aggregatedFunction(matrix, seq, w1, w2, True, [lt])
            aligned = numOfAlignedChars(bitsToStrings(matrix, seq))
            inserted = numOfInsertedIndels(matrix, seq)
            scores.append(score)
            aligns.append(aligned)
            inserts.append(inserted)

            sumScore += score
            sumAligned += aligned
            sumInserted += inserted

            logging.info("\tScore: " + str(score))
            print("\tScore: " + str(score))

            if score > bestScore:
                bestAligned = aligned
                bestPos = pos
                bestScore = score
                bestInserted = inserted
                bestBitString = matrix

            logging.info("\tBest Pos: " + str(bestPos))
            logging.info("\tBest Score: " + str(bestScore))
            print("\tBest Pos: " + str(bestPos))
            print("\tBest Score: " + str(bestScore))

    s = stdev(scores)
    i = stdev(inserts)
    a = stdev(aligns)

    logging.info("Best Score: " + str(bestScore))
    logging.info("Avg Score: " + str(sumScore / 30))
    logging.info("Best Aligned: " + str(bestAligned))
    logging.info("Avg Aligned: " + str(sumAligned / 30))
    logging.info("Best Inserted: " + str(bestInserted))
    logging.info("Avg Inserted:" + str(sumInserted / 30))
    logging.info("St. Dev. Score: " + str(s))
    logging.info("St. Dev. Inserts: " + str(i))
    logging.info("St. Dev. Aligns: " + str(a))
    print("Best Score:", bestScore)
    print("Avg Score:", sumScore / 30)
    print("Best Aligned:", bestAligned)
    print("Avg Aligned:", sumAligned / 30)
    print("Best Inserted:", bestInserted)
    print("Avg Inserted:", sumInserted / 30)
    print("St. Dev. Score: " + str(s))
    print("St. Dev. Inserts: " + str(i))
    print("St. Dev. Aligns: " + str(a))

    for string in bitsToStrings(bestBitString, seq):
        logging.info(string)
        print(string)

    logging.info("Ended " + str(datetime.datetime.now().time()))
    print("Ended " + str(datetime.datetime.now().time()))


logging.basicConfig(filename="ampso " + str(datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S.%f")) + ".txt",
                    level=logging.INFO,
                    format='%(message)s')
'''
test = ["FFABCD", "ABCDFF", "GGABCD", "ABCDGG"]
print("testing")
testAMPSOFuncWeight(test, 0.6, 0.4)
testAMPSOFuncWeight(test, 0.5, 0.5)
testAMPSOFuncWeight(test, 0.3, 0.7)
'''

tests = ["CBCADCAACE", "EACABDCADB", "DABAECBDCD", "DBEACEACCD", "DDABDEEEDE", "EEAECCAAEB", "EABEBCBCCB", "BAADDACDBB"]  # med3
testAMPSOFuncWeight(tests, 0.6, 0.4)
testAMPSOFuncWeight(tests, 0.5, 0.5)
testAMPSOFuncWeight(tests, 0.3, 0.7)

'''
logging.info("test 1")
print("test 1")
testAMPSOFuncWeight(test1, 0.6, 0.4)
testAMPSOFuncWeight(test1, 0.5, 0.5)
testAMPSOFuncWeight(test1, 0.3, 0.7)

logging.info("test 2")
print("test 2")
testAMPSOFuncWeight(test2, 0.6, 0.4)
testAMPSOFuncWeight(test2, 0.5, 0.5)
testAMPSOFuncWeight(test2, 0.3, 0.7)

logging.info("test 3")
print("test 3")
testAMPSOFuncWeight(test3, 0.6, 0.4)
testAMPSOFuncWeight(test3, 0.5, 0.5)
testAMPSOFuncWeight(test3, 0.3, 0.7)

logging.info("test 4")
print("test 4")
testAMPSOFuncWeight(test4, 0.6, 0.4)
testAMPSOFuncWeight(test4, 0.5, 0.5)
testAMPSOFuncWeight(test4, 0.3, 0.7)

logging.info("test 5")
print("test 5")
testAMPSOFuncWeight(test5, 0.6, 0.4)
testAMPSOFuncWeight(test5, 0.5, 0.5)
testAMPSOFuncWeight(test5, 0.3, 0.7)

logging.info("test 6")
print("test 6")
testAMPSOFuncWeight(test6, 0.6, 0.4)
testAMPSOFuncWeight(test6, 0.5, 0.5)
testAMPSOFuncWeight(test6, 0.3, 0.7)

logging.info("test 7")
print("test 7")
testAMPSOFuncWeight(test7, 0.6, 0.4)
testAMPSOFuncWeight(test7, 0.5, 0.5)
testAMPSOFuncWeight(test7, 0.3, 0.7)
'''
