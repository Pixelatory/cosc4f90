import random
import math
import datetime
import logging
import concurrent.futures
import copy
from typing import List
from operator import lt
from util import aggregatedFunction, bitsToStrings, getLongestSeqDict, test1, test2, test3, test4, test5, test6, \
    test7, numOfAlignedChars, numOfInsertedIndels

"""
    BPSO for the MSA Problem
    Nick Aksamit 2020
        
    Acknowledgement goes towards: 
    
    1. A Comparison of Binary Particle Swarm Optimization and Angle Modulated Particle 
    Swarm Optimization for Multiple Sequence Alignment (unknown)
    
    2. Multiple Sequence Alignment Based on a Binary 
    Particle Swarm Optimization Algorithm (Hai-Xia Long, Wen-Bo Xu, Jun Sun, Wen-Juan Ji)
    
    1. Helped me reference that the BPSO works correctly to this person's studies
    2. Helped me reference how the MSA problem should be represented in a BPSO format
"""


def MSABPSO(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, f, w1, w2):
    """The BPSO algorithm fitted for the MSA problem.

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
    :type f: (List[List[int]], List[str], float, float, bool, (int, int) -> bool) -> float
    :param f: fitness function (position vector, sequences, weight coefficient 1, weight coefficient 2)
    :type w1: float
    :param w1: weight coefficient for number of aligned characters
    :type w2: float
    :param w2: weight coefficient for number of leading indels used

    :rtype: List[List[int]]
    :returns: global best position

    Initialization Process:
        Particle positions: each xi = U(0, 1), where U is uniformly distributed random value that's either 0 or 1

        Particle personal best: Same as initialized position

        Particle velocities: each vi = 0

    BPSO:
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
        return f(pos, seq, w1, w2, True, lt)

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
            gBestPos = copy.deepcopy(position)

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
            if fitness(pPositions[i]) > fitness(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = copy.deepcopy(pPositions[i])

            

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(n):
            if fitness(pPositions[i]) > fitness(gBestPos):  # update global best if applicable
                gBestPos = copy.deepcopy(pPositions[i])

        it = it + 1

    return gBestPos


def Sigmoid(x):
    """The classic sigmoid function.

    :type x: float
    :rtype: float
    """
    return 1 / (1 + math.exp(-x))


# ----TESTING AREA----#
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

    logging.info("Started " + str(datetime.datetime.now().time()))
    logging.info("w1: " + str(w1) + " w2: " + str(w2))
    print("Started " + str(datetime.datetime.now().time()))
    print("w1:", w1, "w2:", w2)

    e = []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in range(30):
            e.append(
                executor.submit(MSABPSO, seq, 30, 0.9, 2, 2, 4, 500, float('inf'), 5000, aggregatedFunction, w1, w2))

        for future in concurrent.futures.as_completed(e):
            result = future.result()

            logging.info("A result: " + str(result))
            print("A result: " + str(result))

            score = aggregatedFunction(result, seq, w1, w2, False)
            aligned = numOfAlignedChars(bitsToStrings(result, seq))
            inserted = numOfInsertedIndels(result, seq)

            sumScore += score
            sumAligned += aligned
            sumInserted += inserted

            logging.info("\tScore: " + str(score))
            print("\tScore: " + str(score))

            if score > bestScore:
                bestAligned = aligned
                bestPos = result
                bestScore = score
                bestInserted = inserted

            logging.info("\tBest Pos: " + str(bestPos))
            logging.info("\tBest Score: " + str(bestScore))
            print("\tBest Pos: " + str(bestPos))
            print("\tBest Score: " + str(bestScore))

    logging.info("Best Score: " + str(bestScore))
    logging.info("Avg Score: " + str(sumScore / 30))
    logging.info("Best Aligned: " + str(bestAligned))
    logging.info("Avg Aligned: " + str(sumAligned / 30))
    logging.info("Best Inserted: " + str(bestInserted))
    logging.info("Avg Inserted:" + str(sumInserted / 30))
    print("Best Score:", bestScore)
    print("Avg Score:", sumScore / 30)
    print("Best Aligned:", bestAligned)
    print("Avg Aligned:", sumAligned / 30)
    print("Best Inserted:", bestInserted)
    print("Avg Inserted:", sumInserted / 30)

    for string in bitsToStrings(bestPos, seq):
        logging.info(string)
        print(string)

    logging.info("Ended " + str(datetime.datetime.now().time()))
    print("Ended " + str(datetime.datetime.now().time()))


logging.basicConfig(filename="bpso " + str(datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S.%f")) + ".txt",
                    level=logging.INFO,
                    format='%(message)s')

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
