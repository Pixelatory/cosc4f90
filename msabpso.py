import random
import math
import os
import datetime
import logging
import concurrent.futures
import copy
from typing import List
from operator import eq
from util import aggregatedFunction, bitsToStrings, getLongestSeqDict, test1, test2, test3, test4, test5, test6, test7

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


def MSABPSO(seq, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, f, w1, w2, log):
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
    :type log: bool
    :param log: logging results of the MSABPSO

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

    if log:
        mkdir("MSABPSO logs")
        logging.basicConfig(filename=("MSABPSO logs/msabpso" + str(datetime.datetime.now().time()) + ".log"),
                            format='%(message)s',
                            level=logging.INFO)
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
        :type pos: List[List[int]]
        :param pos: Position vector
        :rtype: float
        :returns: Fitness value of position vector
        """
        return f(pos, seq, w1, w2, True, eq)

    # Just some helper variables to make code more readable
    numOfSeq = len(seq)  # number of sequences
    lSeq = getLongestSeqDict(seq)

    if log:
        logging.info("Number of Sequences: " + str(numOfSeq) + "\n")

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength = math.ceil(lSeq["len"] * 1.2)

    gBestPos = [[0] * colLength] * numOfSeq  # global best position

    if log:
        logging.info("Initial Global Best: " + str(gBestPos))
        logging.info("Initial Global Best Fitness: " + str(fitness(gBestPos)))
        logging.info("\nInitial Particle Values:")

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
            velocity.append([float(0)] * colLength)
            for x in range(colLength - len(seq[j])):
                while True:
                    randNum = random.randint(0, colLength - 1)
                    if position[j][randNum] != 1:
                        position[j][randNum] = 1
                        break

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

    print("w1:", w1, "w2:", w2)

    e = []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in range(30):
            print("Created " + str(i))
            e.append(
                executor.submit(MSABPSO, seq, 30, 0.9, 2, 2, 4, 500, float('inf'), 5000, aggregatedFunction, w1, w2,
                                False))

        for future in concurrent.futures.as_completed(e):
            result = future.result()
            print("A result: " + str(result))
            score = aggregatedFunction(result, seq, w1, w2, False)
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
    for string in bitsToStrings(bestPos, seq):
        print(string)


print("test 1")
testBPSOFuncWeight(test1, 0.6, 0.4)
testBPSOFuncWeight(test1, 0.5, 0.5)
testBPSOFuncWeight(test1, 0.3, 0.7)

print("test 2")
testBPSOFuncWeight(test2, 0.6, 0.4)
testBPSOFuncWeight(test2, 0.5, 0.5)
testBPSOFuncWeight(test2, 0.3, 0.7)

print("test 3")
testBPSOFuncWeight(test3, 0.6, 0.4)
testBPSOFuncWeight(test3, 0.5, 0.5)
testBPSOFuncWeight(test3, 0.3, 0.7)

print("test 4")
testBPSOFuncWeight(test4, 0.6, 0.4)
testBPSOFuncWeight(test4, 0.5, 0.5)
testBPSOFuncWeight(test4, 0.3, 0.7)

print("test 5")
testBPSOFuncWeight(test5, 0.6, 0.4)
testBPSOFuncWeight(test5, 0.5, 0.5)
testBPSOFuncWeight(test5, 0.3, 0.7)

print("test 6")
testBPSOFuncWeight(test6, 0.6, 0.4)
testBPSOFuncWeight(test6, 0.5, 0.5)
testBPSOFuncWeight(test6, 0.3, 0.7)

print("test 7")
testBPSOFuncWeight(test7, 0.6, 0.4)
testBPSOFuncWeight(test7, 0.5, 0.5)
testBPSOFuncWeight(test7, 0.3, 0.7)
