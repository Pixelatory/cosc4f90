import logging
import math
import random
import copy
import concurrent.futures
import os
import datetime
from typing import List, Tuple, Union
from util import getLongestSeqDict, genBitMatrix, numOfAlignedChars, numOfInsertedIndels, bitsToStrings, test1, test2, \
    test3, test4, test5, test6, test7, dominates, theSame, removeCrowdedSolution, updateCrowdingDistances, archiveGuide, \
    addToArchive

"""
    MGPSO for the MSA Problem (using angular modulation)
    Nick Aksamit 2020

    Acknowledgement goes towards:
"""


def MSAMGPSO(seq, genInterval, n, w, c1, c2, c3, k, vmax, vmaxiterlimit, term, maxIter):
    """The MGPSO algorithm fitted for the MSA problem using angular modulation (AM).

        Initialization Process:
            Particle positions: each xi = TODO replace this

            Particle personal best: Same as initialized position

            Particle velocities: each vi = 0

        MGPSO (AM version):
            Topology: Star (gBest)

            PSO Type: Synchronous

            Velocity Update: v(t+1) = w * vi + r1 * c1 * (yi - xi) + l * r2 * c2 * (ŷi - xi) + (1 - l) * r3 * c3 * (
            TODO,
            where r1, r2, r3 are uniformly random values between [0,1]

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
            :type term: List[float]
            :param term: termination criteria for each function (set to float('inf') for no fitness termination)
            :type maxIter: int
            :param maxIter: maximum iteration limit (> 0)
            :rtype: List[List[List[float], List[List[int]], float]]

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
    f = [numOfAlignedChars, numOfInsertedIndels]  # objective functions
    pPositions: List[List[List[float]]] = []  # particle positions
    pPersonalBests: List[List[List[float]]] = []  # particle personal best position
    pVelocities: List[List[List[float]]] = []  # particle velocities
    pBitStrings: List[List[List[List[int]]]] = []  # particle bit strings
    sArchive: List[List[List[float], List[List[int]], float]] = []  # swarm archive (the pareto front)

    lSeq = getLongestSeqDict(seq)  # Longest sequence value dictionary

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength: int = math.ceil(lSeq["len"] * 1.2)

    # Initializing gBest
    # Each subswarm has its own gBest pos and bitmatrix
    gBest = []

    # Initializing the lambda parameter
    l: float = 0

    for i in range(len(f)):
        gBest.append({})
        gBest[i]["pos"] = [0, 0, 0, 0]
        gBest[i]["bitstring"] = genBitMatrix(gBest[i]["pos"], seq, colLength, genInterval)

    # Initializing the particles of each subswarm
    for i in range(len(f)):
        pPositions.append([])
        pPersonalBests.append([])
        pVelocities.append([])
        pBitStrings.append([])

        for j in range(n):
            position: List[float] = []
            velocity: List[float] = [0] * 4

            position.append(random.uniform(-1, 1))  # coefficient a
            position.append(random.uniform(0, 1))  # coefficient b
            position.append(random.uniform(0, 1))  # coefficient c
            position.append(random.uniform(-0.7, -0.9))  # coefficient d

            bitstring = genBitMatrix(position, seq, colLength, genInterval)

            pPositions[i].append(position)
            pPersonalBests[i].append(position)
            pVelocities[i].append(velocity)
            pBitStrings[i].append(bitstring)

            if i == 0:  # numOfAlignedChars
                strings = bitsToStrings(bitstring, seq)
                gBestStrings = bitsToStrings(gBest[i]["bitstring"], seq)
                if f[i](strings) > f[i](gBestStrings):
                    gBest[i]["pos"] = copy.deepcopy(position)
                    gBest[i]["bitstring"] = copy.deepcopy(bitstring)
            else:  # numOfInsertedIndels
                if f[i](bitstring, seq) < f[i](gBest[i]["bitstring"], seq):
                    gBest[i]["pos"] = copy.deepcopy(position)
                    gBest[i]["bitstring"] = copy.deepcopy(bitstring)

    def multiThreaded(i, j):
        nonlocal pPositions, pVelocities, pBitStrings, pPersonalBests, gBest, vmaxiterlimit, vmax, seq, colLength, \
            genInterval, f, k, l
        # Within each subswarm, update each particle's velocity, position, and personal best
        # r1 and r2 are ~ U (0,1)
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        a = archiveGuide(seq, sArchive, 1, 2, k)

        # update velocity and positions in every dimension
        for k in range(4):
            pVelocities[i][j][k] = w * pVelocities[i][j][k] + \
                                   r1 * c1 * (pPersonalBests[i][j][k] - pPositions[i][j][k]) + \
                                   l * r2 * c2 * (gBest[i]["pos"][k] - pPositions[i][j][k]) + \
                                   (1 - l) * r3 * c3 * (a[k] - pPositions[i][j][k])

            # velocity clamping
            if vmaxiterlimit < it:
                if pVelocities[i][j][k] > vmax:
                    pVelocities[i][j][k] = vmax
                elif pVelocities[i][j][k] < -vmax:
                    pVelocities[i][j][k] = -vmax

            pPositions[i][j][k] += pVelocities[i][j][k]

        bitstring = genBitMatrix(pPositions[i][j], seq, colLength, genInterval)

        # update personal best if applicable
        if i == 0:  # numOfAlignedChars
            strings = bitsToStrings(bitstring, seq)
            pBestStrings = bitsToStrings(pBitStrings[i][j], seq)
            if f[i](strings) > f[i](pBestStrings):
                pPersonalBests[i][j] = copy.deepcopy(pPositions[i][j])
                pBitStrings[i][j] = copy.deepcopy(bitstring)
        else:  # numOfInsertedIndels
            if f[i](bitstring, seq) < f[i](pBitStrings[i][j], seq):
                pPersonalBests[i][j] = copy.deepcopy(pPositions[i][j])
                pBitStrings[i][j] = copy.deepcopy(bitstring)

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter:
        # print(str(it) + " " + str(it / maxIter) + "%")

        # if the termination criteria is met, then stop the PSO and return values
        #  numOfAlignedChars                                            numOfInsertedIndels
        if f[0](bitsToStrings(gBest[0]["bitstring"], seq)) > term[0] or f[1](gBest[1]["bitstring"], seq) < term[1]:
            return sArchive

        for i in range(len(f)):
            for j in range(n):
                sArchive = addToArchive(seq, sArchive, (pPositions[i][j], pBitStrings[i][j]), 1, 2, len(f) * n)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            for i in range(len(f)):
                for j in range(n):
                    executor.submit(multiThreaded, i, j)

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(len(f)):
            for j in range(n):
                # update global best if applicable
                if i == 0:  # numOfAlignedChars
                    strings = bitsToStrings(pBitStrings[i][j], seq)
                    gBestStrings = bitsToStrings(gBest[i]["bitstring"], seq)
                    if f[i](strings) > f[i](gBestStrings):
                        gBest[i]["pos"] = copy.deepcopy(pPositions[i][j])
                        gBest[i]["bitstring"] = copy.deepcopy(pBitStrings[i][j])
                else:  # numOfInsertedIndels
                    if f[i](pBitStrings[i][j], seq) < f[i](gBest[i]["bitstring"], seq):
                        gBest[i]["pos"] = copy.deepcopy(pPositions[i][j])
                        gBest[i]["bitstring"] = copy.deepcopy(pBitStrings[i][j])

        # update the lambda parameter (linearly increasing)
        l += 1 / maxIter

        it = it + 1

    return sArchive


def testing(seqs, i, iterations):
    print("\nTest " + str(i))
    logging.info("\nTest " + str(i))
    logging.info(seqs)
    logging.info("genInterval = [-2.0, 2.0]")
    logging.info("n = 30")
    logging.info("w = 0.729844")
    logging.info("c1 = c2 = c3 = 1.49618")
    logging.info("k = 3")
    logging.info("vmax = infinite")
    logging.info("term = maxIter")
    logging.info("maxIter = " + str(iterations))

    print(str(iterations) + " iterations:")
    t = MSAMGPSO(seqs, [-2.0, 2.0], 30, 0.729844, 1.49618, 1.49618, 1.49618, 3, float('inf'), 500,
                 [float('inf'), -float('inf')], iterations)
    for res in t:
        logging.info(res[1])
        print(res[1])
        logging.info(res[0])
        print(res[0])
        for string in bitsToStrings(res[1], seqs):
            logging.info(string)
            print(string)
        logging.info(numOfAlignedChars(bitsToStrings(res[1], seqs)))
        print(numOfAlignedChars(bitsToStrings(res[1], seqs)))
        logging.info(numOfInsertedIndels(res[1], seqs))
        print(numOfInsertedIndels(res[1], seqs))


logging.basicConfig(filename="mpampso " + str(datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S.%f")) + ".txt",
                    level=logging.INFO,
                    format='%(message)s')

AB000177 = "gaccatatgattgacgcctatgtcaatctctacactacattgctggaaagcaaatcctgagagatgctacccccgccgttgctgcgggggccaacgcgttaatgccgattcttcagattatcaatcacttctccgagatccagcccctgatcctgcaacagcaccagcaggtgatacaccaaatcagatgcctcattcttcagctcaaagcggtcatttaccgttgcggccagtgcggttt"
AB000178 = "gaccatatgattgacgcctatgtcaatctctacactacattgctggaaagcaaatcctgagagatgctacccccgccgttgctgcgggggccaatgcgttaatgccgattcttcagattatcaatcacttctccgagatccagcccctgatcctgtaacagcaccagcaggtgatacatcaaatcagatgcctcgttggtcagctcaaagcggtcatgtaccgttggtgccagtgcggttt"
AB000179 = "gaccatatgattgacgcctatgtcaatctctacactacattgctggaaagcaaatcctgagagatgctacccccgccgttgctgcgggggccaatgcgttaatgccgattcttcagattatcaatcacttctccgagatccagcccctgatcctgtaacagcaccagcaggtgatacatcaaatcagatgcctcgttggtcagctcaaagcggtcatgtaccgtcgcggccagtgcagttt"
AB000180 = "gaccatatgattgacgcctatgtcaatctctacactacattgctggaaagcaaatcctgagagatgctacccccgccgttgctgcgggggccaacgcgttaatgccgattcttcagattatcaatcacttctccgagatccagcccttgatcctgcaacaacaccagcaggtgatacatcaaatcagatgcctcgttggtcagttcaaagcggtcatgtactgtcgctgccagtgcggttt"

strs = [AB000177, AB000178, AB000179, AB000180]

testing(strs, 1, 1000)
testing(strs, 2, 2500)
testing(strs, 3, 5000)
testing(strs, 3, 7500)
testing(strs, 4, 10000)

'''
testing(test1, 1)
testing(test2, 2)
testing(test3, 3)
testing(test4, 4)
testing(test5, 5)
testing(test6, 6)
testing(test7, 7)

#os.system('shutdown -s') # shutdown computer'''
