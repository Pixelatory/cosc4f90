import math
import random
import copy
from typing import List, Tuple
from util import getLongestSeqDict, genBitMatrix, numOfAlignedChars, numOfInsertedIndels, bitsToStrings

"""
    MGPSO for the MSA Problem
    Nick Aksamit 2020

    Acknowledgement goes towards:
"""


# TODO: find a way to generalize the objective functions if possible. If not then I mean hard-coded is still alright.

def MSAMGPSO(seq, genInterval, coefLimit, n, w, c1, c2, c3, l, k, vmax, vmaxiterlimit, term, maxIter):
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

    # Initialize the data containers
    f = [numOfAlignedChars, numOfInsertedIndels]  # objective functions
    pPositions: List[List[List[float]]] = []  # particle positions
    pPersonalBests: List[List[List[float]]] = []  # particle personal best position
    pVelocities: List[List[List[float]]] = []  # particle velocities
    pBitStrings: List[List[List[List[int]]]] = []  # particle bit strings
    sArchive: List[List[List[float], List[List[int], float]]] = []  # swarm archive (the pareto front)

    def dominates(bm1, bm2):
        better = False
        for x in range(len(f)):
            if x == 0:  # numOfAlignedChars
                strings1 = bitsToStrings(bm1, seq)
                strings2 = bitsToStrings(bm2, seq)
                res1 = numOfAlignedChars(strings1)
                res2 = numOfAlignedChars(strings2)
                if res1 < res2:
                    return False
                elif res1 > res2:
                    better = True
            else:  # numOfInsertedIndels
                res1 = numOfInsertedIndels(bm1, seq)
                res2 = numOfInsertedIndels(bm2, seq)
                if res1 > res2:
                    return False
                elif res1 < res2:
                    better = True
        return better

    def addToArchive(x):
        """If possible, adds solution x to archive a.

        :type x: Tuple[List[float], List[List[int]]]
        :rtype: None
        """

        global sArchive
        aDominated = []  # all the archive components that are dominated by x

        # First, check that this new solution dominates every archive solution
        for s in sArchive:
            if not dominates(x[1], s[1]):
                return
            else:
                aDominated.append(s)

        # When the function reaches here, it's safe to say the solution dominates.
        # So, let's add it to archive.
        sArchive.append(copy.deepcopy([x[0], x[1], 0.0]))

        # Next, after adding it remove all the archive elements that x dominates
        sArchive = [x for x in sArchive if x not in aDominated]

        if len(sArchive) > len(f) * n:
            removeCrowdedSolution()

    def updateCrowdingDistance():
        """Calculates the crowding distance of each archive solution."""
        for i in range(len(f)):
            if i == 0:  # numOfAlignedChars
                sArchive.sort(key=lambda x: numOfAlignedChars(bitsToStrings(x[1], seq)), reverse=True)
            else:  # numOfInsertedIndels
                sArchive.sort(key=lambda x: numOfInsertedIndels(x[1], seq))

            sArchive[0][2] = sArchive[len(sArchive) - 1][2] = float('inf')  # set first and last to infinite distance

            for j in range(1, len(sArchive) - 2):
                if sArchive[j][2] != float('inf'):
                    sArchive[j][2] += sArchive[j+1][2] - sArchive[j-1][2]

    def removeCrowdedSolution():
        """Removes the most crowded solution in the archive."""
        updateCrowdingDistance()
        minIdx = 0

        for i in range(len(sArchive)):
            if sArchive[i][2] < sArchive[minIdx][2]:
                minIdx = i

        del sArchive[minIdx]

    def archiveGuide():
        """Uses tournament selection where k is the number of particles to choose.

        Out of the k possible particles randomly selcted, the least crowded particle wins the tournament.
        """



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

    # This is where the iterations begin
    it = 0  # iteration count
    while it < maxIter:
        # if the termination criteria is met, then stop the PSO and return values
        if f[0](bitsToStrings(gBest[0]["bitstring"], seq)) > term[0]:
            return True
        elif f[1](gBest[1]["bitstring"], seq) > term[1]:
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

        it = it + 1
