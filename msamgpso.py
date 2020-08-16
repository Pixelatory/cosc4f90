import math
import random
import copy
from typing import List, Tuple
from operator import lt
from util import getLongestSeqDict, genBitMatrix

"""
    MGPSO for the MSA Problem
    Nick Aksamit 2020

    Acknowledgement goes towards:
"""


def MSAMGPSO(seq, genInterval, coefLimit, n, w, c1, c2, c3, l, vmax, vmaxiterlimit, term, maxIter, f, w1, w2):
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
    elif type(f) is not list:
        raise Exception("f must be a list of callable functions")
    elif len(term) != len(f):
        raise Exception("termination list must be same size as function list")

    # sort coefLimit and genInterval, just makes life easier
    if coefLimit[0] > coefLimit[1]:
        tmp = coefLimit[0]
        coefLimit[0] = coefLimit[1]
        coefLimit[1] = tmp

    if genInterval[0] > genInterval[1]:
        tmp = genInterval[0]
        genInterval[0] = genInterval[1]
        genInterval[1] = tmp

    # Initialize the data containers
    pPositions: List[List[List[float]]] = []  # particle positions
    pPersonalBests: List[List[List[float]]] = []  # particle personal best position
    pVelocities: List[List[List[float]]] = []  # particle velocities
    pBitStrings: List[List[List[List[int]]]] = []  # particle bit strings

    # swarm archive
    sArchive = []  # TODO: figure out the archive

    def addToArchive(x):
        """If possible, adds solution x to archive a.

        :type x: Tuple[List[float], List[List[int]]]
        """
        for s in sArchive:


    def fitness(bitmatrix, f):
        """
        To test fitness in the AMPSO, first you use the position vector as the coefficients
        of the angular modulation formula. Then, sample random values within genInterval with
        the coefficients and use these values with the gen function. If the gen function
        returns a value > 0, the bit is 1, otherwise 0.

        :type bitmatrix: List[List[int]]
        :param bitmatrix: Two-dimensional binary matrix
        :type f: (List[List[int]], List[str], float, float, bool, (int, int) -> bool) -> float
        :rtype: float
        :returns: Fitness value of bit string
        """

        return f(bitmatrix, seq, w1, w2, True, lt)

    lSeq = getLongestSeqDict(seq)  # Longest sequence value dictionary

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength: int = math.ceil(lSeq["len"] * 1.2)

    # Initializing gBest
    # Each subswarm has its own gBest pos and bitmatrix
    gBest = []

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

            if fitness(bitstring, f[i]) > fitness(gBest[i]["bitstring"], f[i]):
                gBest[i]["pos"] = copy.deepcopy(position)
                gBest[i]["bitstring"] = copy.deepcopy(bitstring)

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter:
        # if the termination criteria is met, then stop the PSO and return values
        for i in range(len(gBest)):
            if fitness(gBest[i]["bitstring"], f[i]) > term[i]:
                return True

        # Update each particle's velocity, position, and personal best
        for i in range(len(f)):
            for j in range(n):
                # r1 and r2 are ~ U (0,1)
                r1 = random.random()
                r2 = random.random()
                r3 = random.random()

                # update velocity and positions in every dimension
                for k in range(4):
                    pVelocities[i][j][k] = w * pVelocities[i][j][k] + \
                                           r1 * c1 * (pPersonalBests[i][j][k] - pPositions[i][j][k]) + \
                                           l * r2 * c2 * (gBest[i]["pos"][k] - pPositions[i][j][k]) + \
                                           (1 - l) * r3 * c3 * (sArchive[i][j][k] - pPositions[i][j][k])

                    # velocity clamping
                    if vmaxiterlimit < it:
                        if pVelocities[i][j][k] > vmax:
                            pVelocities[i][j][k] = vmax
                        elif pVelocities[i][j][k] < -vmax:
                            pVelocities[i][j][k] = -vmax

                    pPositions[i][j][k] = pPositions[i][j][k] + pVelocities[i][j][k]

                    '''if pPositions[i][j] > coefLimit[1]:
                        pPositions[i][j] = coefLimit[1]
                    elif pPositions[i][j] < coefLimit[0]:
                        pPositions[i][j] = coefLimit[0]'''

                bitstring = genBitMatrix(pPositions[i][j], seq, colLength, genInterval)

                # update personal best if applicable
                if fitness(bitstring, f[i]) > fitness(pBitStrings[i][j], f[i]):  # update personal best if applicable
                    pPersonalBests[i][j] = copy.deepcopy(pPositions[i][j])
                    pBitStrings[i][j] = copy.deepcopy(bitstring)

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(len(f)):
            for j in range(n):
                # update global best if applicable
                if fitness(pBitStrings[i][j], f[i]) > fitness(gBest[i]["bitstring"], f[i]):
                    gBest[i]["pos"] = copy.deepcopy(pPositions[i][j])
                    gBest[i]["bitstring"] = copy.deepcopy(pBitStrings[i][j])

        it = it + 1
