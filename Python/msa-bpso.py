import random
import math
import datetime
import logging
import concurrent.futures

from copy import deepcopy
from statistics import stdev
from typing import List
from operator import lt, gt
from util import aggregatedFunction, bitsToStrings, getLongestSeqDict, test1, test2, test3, test4, test5, test6, \
    test7, numOfAlignedChars, numOfInsertedIndels, Sigmoid

"""
    pso.BPSO for the MSA Problem
    Nick Aksamit 2020
        
    Acknowledgement goes towards: 
    
    1. Multiple Sequence Alignment Based on a Binary 
    Particle Swarm Optimization Algorithm (Hai-Xia Long, Wen-Bo Xu, Jun Sun, Wen-Juan Ji)
    
    1. Helped me reference how the MSA problem should be represented in a pso.BPSO format
    
    Note of the imported functions to reduce confusion in code:
    - deepcopy (from copy)
    - stdev (from stdev)
    - List (from typing)
    - lt (from operator)
    - aggregatedFunction (from util)
    - bitsToString (from util)
    - getLongestSeqDict (from util)
    - test1 (from util)
    - test2 (from util)
    - test3 (from util)
    - test4 (from util)
    - test5 (from util)
    - test6 (from util)
    - test7 (from util)
    - numOfAlignedChars (from util)
    - numOfInsertedIndels (from util)
    - Sigmoid (from util)
"""


def MSABPSO(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, f, w1, w2, ops):
    """The pso.BPSO algorithm fitted for the MSA problem.

    :type seq: List[str]
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
    :type vmaxiterlimit: int
    :param vmaxiterlimit: Maximum iteration limit of clamping velocity values
    :type term: float
    :param term: termination criteria (set to float('inf') for no fitness termination)
    :type maxIter: int
    :param maxIter: maximum iteration limit (> 0)
    :type f: (List[List[int]], List[str], float, float, bool, List[(int, int) -> bool]) -> float
    :param f: fitness function (position vector, sequences, weight coefficient 1, weight coefficient 2)
    :type w1: float
    :param w1: weight coefficient for number of aligned characters
    :type w2: float
    :param w2: weight coefficient for number of leading indels used
    :param ops: operators that will check for infeasibility (See util.py -> infeasible)
    :type ops: List[(float, float) -> bool]

    :rtype: (List[List[int]], int)
    :returns: (global best position, numOfInfeasibleSols)

    Initialization Process:
        Particle positions: each xi = U(0, 1), where U is uniformly distributed random value that's either 0 or 1

        Particle personal best: Same as initialized position

        Particle velocities: each vi = 0

    pso.BPSO:
        Topology: Star (gBest)

        PSO Type: Synchronous

        Velocity Update: v(t+1) = w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi),
        where r1 and r2 are uniformly random values between [0,1]

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
    pPositions: List[List[List[int]]] = []
    pPersonalBests: List[List[List[int]]] = []
    pVelocities: List[List[List[float]]] = []
    numOfInfeasibleSols: int = 0

    def fitness(pos):
        """
        :type pos: List[List[int]]
        :param pos: Position vector
        :rtype: float
        :returns: Fitness value of position vector
        """
        return f(pos, seq, w1, w2, True, ops)

    # Just some helper variables to make code more readable
    numOfSeq = len(seq)  # number of sequences
    lSeq = getLongestSeqDict(seq)

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength = math.ceil(lSeq["len"] * 1.2)

    gBestPos = [[0] * colLength] * numOfSeq  # global best position

    # Initializing the particles of swarm
    for i in range(n):
        position: List[List[int]] = []
        velocity: List[List[float]] = [[float(0)] * colLength] * len(seq)

        # This initiates the position and velocity, which are
        # both a matrix of size [numOfSequences x colLength].
        # Additionally it randomly puts indels in a random spot,
        # but only just the right amount of them so the initial
        # position is feasible.
        for j in range(numOfSeq):
            position.append([0] * colLength)
            for x in range(colLength - len(seq[j])):
                while True:
                    randNum = random.randint(0, colLength - 1)
                    if position[j][randNum] != 1:
                        position[j][randNum] = 1
                        break

        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)

        if fitness(position) > fitness(gBestPos):
            gBestPos = deepcopy(position)

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter and fitness(gBestPos) < term:

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
                                           r2 * c2 * (gBestPos[j][x] - pPositions[i][j][x])

                    # velocity clamping
                    if vmaxiterlimit < it:
                        if pVelocities[i][j][x] > vmax:
                            pVelocities[i][j][x] = vmax
                        elif pVelocities[i][j][x] < -vmax:
                            pVelocities[i][j][x] = -vmax

                    probability = Sigmoid(pVelocities[i][j][x])
                    pPositions[i][j][x] = 1 if random.uniform(0, 1) < probability else 0

            # update personal best if applicable
            tmpScore = fitness(pPositions[i])

            # Infeasible solution, increment counter and don't attempt adding to personal best
            if tmpScore == -float('inf'):
                numOfInfeasibleSols += 1
                continue

            if tmpScore > fitness(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = deepcopy(pPositions[i])

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(n):
            pFit = fitness(pPositions[i])

            # Infeasible solution, don't attempt adding to global best
            if pFit == -float('inf'):
                continue

            if pFit > fitness(gBestPos):  # update global best if applicable
                gBestPos = deepcopy(pPositions[i])

        it = it + 1

    return gBestPos, numOfInfeasibleSols


# ----TESTING AREA----#
def testBPSOFuncWeight(seq, w1, w2):
    """Just a testing function for the pso.BPSO on an MSA problem  (dynamic w1 and w2 values)\n
    n = 30\n
    w = 0.9\n
    c1 = c2 = 2\n
    vmax = 2\n
    vmaxiterlimit = 500\n
    term = float('inf')\n
    iter = 5000\n
    f = aggregatedFunction\n
    log = False\n

    :type seq: List[str]
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
                executor.submit(MSABPSO, seq, 30, 0.99, 2, 2, 4, 500, float('inf'), 5000, aggregatedFunction, w1, w2, [lt]))

        for future in concurrent.futures.as_completed(e):
            result = future.result()
            matrix = result[0]
            numOfInfeasibleSols = result[1]

            # Just logging and printing to screen
            logging.info("A result: " + str(bestPos))
            logging.info("Infeasible Sol: " + str((numOfInfeasibleSols / (30 * 5000)) * 100) + "%")
            print("A result: " + str(matrix))
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
                bestPos = matrix
                bestScore = score
                bestInserted = inserted

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
    print("Standard Deviation: " + str(s))
    print("St. Dev. Score: " + str(s))
    print("St. Dev. Inserts: " + str(i))
    print("St. Dev. Aligns: " + str(a))

    for string in bitsToStrings(bestPos, seq):
        logging.info(string)
        print(string)

    logging.info("Ended " + str(datetime.datetime.now().time()))
    print("Ended " + str(datetime.datetime.now().time()))


logging.basicConfig(filename="bpso " + str(datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S.%f")) + ".txt",
                    level=logging.INFO,
                    format='%(message)s')
'''
test = ["FFABCD", "ABCDFF", "GGABCD", "ABCDGG"]  # basic1
print("testing")
testBPSOFuncWeight(test, 0.6, 0.4)
testBPSOFuncWeight(test, 0.5, 0.5)
testBPSOFuncWeight(test, 0.3, 0.7)
'''

tests = ["CBCADCAACE", "EACABDCADB", "DABAECBDCD", "DBEACEACCD", "DDABDEEEDE", "EEAECCAAEB", "EABEBCBCCB", "BAADDACDBB"]  # med3
testBPSOFuncWeight(tests, 0.6, 0.4)
testBPSOFuncWeight(tests, 0.5, 0.5)
testBPSOFuncWeight(tests, 0.3, 0.7)

'''
logging.info("test 1")
print("test 1")
testBPSOFuncWeight(test1, 0.6, 0.4)
testBPSOFuncWeight(test1, 0.5, 0.5)
testBPSOFuncWeight(test1, 0.3, 0.7)

logging.info("test 2")
print("test 2")
testBPSOFuncWeight(test2, 0.6, 0.4)
testBPSOFuncWeight(test2, 0.5, 0.5)
testBPSOFuncWeight(test2, 0.3, 0.7)

logging.info("test 3")
print("test 3")
testBPSOFuncWeight(test3, 0.6, 0.4)
testBPSOFuncWeight(test3, 0.5, 0.5)
testBPSOFuncWeight(test3, 0.3, 0.7)

logging.info("test 4")
print("test 4")
testBPSOFuncWeight(test4, 0.6, 0.4)
testBPSOFuncWeight(test4, 0.5, 0.5)
testBPSOFuncWeight(test4, 0.3, 0.7)

logging.info("test 5")
print("test 5")
testBPSOFuncWeight(test5, 0.6, 0.4)
testBPSOFuncWeight(test5, 0.5, 0.5)
testBPSOFuncWeight(test5, 0.3, 0.7)

logging.info("test 6")
print("test 6")
testBPSOFuncWeight(test6, 0.6, 0.4)
testBPSOFuncWeight(test6, 0.5, 0.5)
testBPSOFuncWeight(test6, 0.3, 0.7)

logging.info("test 7")
print("test 7")
testBPSOFuncWeight(test7, 0.6, 0.4)
testBPSOFuncWeight(test7, 0.5, 0.5)
testBPSOFuncWeight(test7, 0.3, 0.7)
'''