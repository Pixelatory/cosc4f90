import concurrent.futures
import copy
import math
import random
from typing import List
from util import numOfAlignedChars, numOfInsertedIndels, getLongestSeqDict, genBitMatrix, bitsToStrings, archiveGuide, \
    addToArchive

"""
    MGPSO for the MSA Problem (using binary representation)
    Nick Aksamit 2020

    Acknowledgement goes towards:
"""


def MSAMGPSO(seq, n, w, c1, c2, c3, l, k, vmax, vmaxiterlimit, term, maxIter):
    # Checking for trivial errors first
    if n < 1:
        raise Exception("Swarm size cannot be < 1")
    elif len(seq) < 2:
        raise Exception("Number of sequences cannot be < 2")
    elif maxIter < 1:
        raise Exception("maxIter cannot be < 1")
    elif maxIter == float('inf') and term == float('inf'):
        raise Exception("Maximum iterations and termination fitness are both infinite!")

    # Initialize the data containers
    f = [numOfAlignedChars, numOfInsertedIndels]  # objective functions
    pPositions: List[List[List[List[int]]]] = []  # particle positions
    pPersonalBests: List[List[List[List[int]]]] = []  # particle personal best position
    pVelocities: List[List[List[List[float]]]] = []  # particle velocities
    sArchive: List[List[List[List[int], float]]] = []  # swarm archive (the pareto front)

    lSeq = getLongestSeqDict(seq)  # Longest sequence value dictionary

    # The position column length is 20% greater than the total length of longest sequence (rounded up)
    colLength: int = math.ceil(lSeq["len"] * 1.2)

    # Initializing gBest
    # Each subswarm has its own gBest position
    gBest = []

    for i in range(len(f)):
        gBest.append([[0] * colLength] * len(seq))

    # Initializing the particles of each subswarm
    for i in range(len(f)):
        newSwarmPositions = []
        newSwarmPersonalBests = []
        newSwarmVelocities = []

        for j in range(n):
            position: List[List[int]] = []
            velocity: List[List[float]] = [[float(0)] * colLength] * len(seq)

            # This initiates the position and velocity, which are
            # both a matrix of size [numOfSequences x colLength].
            # Additionally it randomly puts indels in a random spot,
            # but only just the right amount of them so the initial
            # position is feasible.
            for x in range(len(seq)):
                position.append([0] * colLength)
                for y in range(colLength - len(seq[x])):
                    while True:
                        randNum = random.randint(0, colLength - 1)
                        if position[x][randNum] != 1:
                            position[x][randNum] = 1
                            break

            newSwarmPositions.append(position)
            newSwarmPersonalBests.append(position)
            newSwarmVelocities.append(velocity)

            if i == 0:  # numOfAlignedChars
                strings = bitsToStrings(position, seq)
                gBestStrings = bitsToStrings(gBest[i], seq)
                if f[i](strings) > f[i](gBestStrings):
                    gBest[i] = copy.deepcopy(position)
            else:  # numOfInsertedIndels
                if f[i](position, seq) < f[i](gBest[i], seq):
                    gBest[i] = copy.deepcopy(position)

    def multiThreaded(i, j):
        nonlocal pPositions, pVelocities, pPersonalBests, gBest, vmaxiterlimit, vmax, seq, colLength, f, k
        # Within each subswarm, update each particle's velocity, position, and personal best
        # r1 and r2 are ~ U (0,1)
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        a = archiveGuide(seq, sArchive, 1, 2, k)

        # update velocity and positions in every dimension
        for x in range(len(seq)):
            for y in range(len(seq[x])):
                pVelocities[i][j][x][y] = w * pVelocities[i][j][x][y] + \
                                          r1 * c1 * (pPersonalBests[i][j][x][y] - pPositions[i][j][x][y]) + \
                                          l * r2 * c2 * (gBest[i][x][y] - pPositions[i][j][x][y]) + \
                                          (1 - l) * r3 * c3 * (a[x][y] - pPositions[i][j][x][y])

                # velocity clamping
                if vmaxiterlimit < it and (pVelocities[i][j][x][y] > vmax or pVelocities[i][j][x][y] < vmax):
                    pVelocities[i][j][x][y] = vmax

                pPositions[i][j][x][y] += pVelocities[i][j][x][y]

        # update personal best if applicable
        if i == 0:  # numOfAlignedChars
            strings = bitsToStrings(pPositions[i][j], seq)
            pBestStrings = bitsToStrings(pPersonalBests[i][j], seq)
            if f[i](strings) > f[i](pBestStrings):
                pPersonalBests[i][j] = copy.deepcopy(pPositions[i][j])
        else:  # numOfInsertedIndels
            if f[i](pPositions[i][j], seq) < f[i](pPersonalBests[i][j], seq):
                pPersonalBests[i][j] = copy.deepcopy(pPositions[i][j])

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter:
        # if the termination criteria is met, then stop the PSO and return values
        #  numOfAlignedChars                               numOfInsertedIndels
        if f[0](bitsToStrings(gBest[0], seq)) > term[0] or f[1](gBest[1], seq) < term[1]:
            return sArchive

        for i in range(len(f)):
            for j in range(n):
                sArchive = addToArchive(seq, sArchive, pPositions[i][j], 0, 1, len(f) * n)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            for i in range(len(f)):
                for j in range(n):
                    executor.submit(multiThreaded, i, j)

        # update the global best after all positions were changed (synchronous PSO)
        for i in range(len(f)):
            for j in range(n):
                # update global best if applicable
                if i == 0:  # numOfAlignedChars
                    strings = bitsToStrings(pPositions[i][j], seq)
                    gBestStrings = bitsToStrings(gBest[i], seq)
                    if f[i](strings) > f[i](gBestStrings):
                        gBest[i] = copy.deepcopy(pPositions[i][j])
                else:  # numOfInsertedIndels
                    if f[i](pPositions[i][j], seq) < f[i](gBest[i], seq):
                        gBest[i] = copy.deepcopy(pPositions[i][j])

        it = it + 1

    return sArchive
