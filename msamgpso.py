import math
import random
import copy
import concurrent.futures
from typing import List, Tuple
from util import getLongestSeqDict, genBitMatrix, numOfAlignedChars, numOfInsertedIndels, bitsToStrings, test1, test2, \
    test3, test4, test5, test6, test7

"""
    MGPSO for the MSA Problem
    Nick Aksamit 2020

    Acknowledgement goes towards:
"""


def MSAMGPSO(seq, genInterval, n, w, c1, c2, c3, l, k, vmax, vmaxiterlimit, term, maxIter):
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

    def theSame(x, y):
        """Checks if two sets of position and bitstring values are the same.

        Assumes that the position vector and bitstrings are of the same length.

        :type x: Tuple[List[float], List[List[int]]]
        :type y: Tuple[List[float], List[List[int]]]
        :rtype: bool
        """

        for i in range(len(x[0])):
            if x[0][i] != y[0][i]:
                return False

        for i in range(len(x[1])):
            for j in range(len(x[1][i])):
                if x[1][i][j] != y[1][i][j]:
                    return False

        return True

    def addToArchive(x):
        """Adds solution x to archive a if x is not dominated by any archive solutions.

        After adding, if any particles in the archive are now dominated, then they are removed.

        Additionally, if the archive is full then the most crowded solution is removed.

        :type x: Tuple[List[float], List[List[int]]]
        :rtype: None
        """
        nonlocal sArchive
        aDominated = []  # all the archive components that are dominated by x

        # First, check that this new solution dominates every archive solution,
        # and that the values (bitstring and position) aren't repeated
        for s in sArchive:
            if dominates(s[1], x[1]):
                return
            elif dominates(x[1], s[1]):
                aDominated.append(s)

            if theSame((s[0], s[1]), x):
                return

        # When the function reaches here, it's safe to say the solution dominates.
        # So, let's add it to archive.
        sArchive.append(copy.deepcopy([x[0], x[1], 0.0]))

        # Next, after adding it remove all the archive elements that x dominates
        sArchive = [x for x in sArchive if x not in aDominated]

        if len(sArchive) > len(f) * n:
            removeCrowdedSolution()

    def updateCrowdingDistances():
        """Calculates the crowding distance of each archive solution."""
        for i in range(len(f)):
            if i == 0:  # numOfAlignedChars
                sArchive.sort(key=lambda x: -numOfAlignedChars(bitsToStrings(x[1], seq)))
            else:  # numOfInsertedIndels
                sArchive.sort(key=lambda x: numOfInsertedIndels(x[1], seq))

            sArchive[0][2] = sArchive[len(sArchive) - 1][2] = float('inf')  # set first and last to infinite distance

            for j in range(1, len(sArchive) - 2):
                if sArchive[j][2] != float('inf'):
                    sArchive[j][2] += sArchive[j + 1][2] - sArchive[j - 1][2]

    def removeCrowdedSolution():
        """Removes the most crowded solution in the archive."""
        updateCrowdingDistances()
        minIdx = 0

        for i in range(len(sArchive)):
            if sArchive[i][2] < sArchive[minIdx][2]:
                minIdx = i

        del sArchive[minIdx]

    def archiveGuide():
        """Uses tournament selection where k is the number of particles to choose.

        Out of the k possible particles randomly selected, the least crowded particle wins the tournament.

        :rtype: List[float]
        :returns: Position vector
        """
        updateCrowdingDistances()

        idx = random.randint(0, len(sArchive) - 1)

        for i in range(k - 1):
            tmp = random.randint(0, len(sArchive) - 1)
            if sArchive[tmp][2] > sArchive[idx][2]:
                idx = tmp

        return sArchive[idx][0]

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

    def multiThreaded(i, j):
        nonlocal pVelocities, pBitStrings, pPersonalBests, gBest, vmaxiterlimit, vmax, seq, colLength, genInterval, f
        # Within each subswarm, update each particle's velocity, position, and personal best
        # r1 and r2 are ~ U (0,1)
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        a = archiveGuide()

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

            pPositions[i][j][k] = pPositions[i][j][k] + pVelocities[i][j][k]

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
            return True

        for i in range(len(f)):
            for j in range(n):
                addToArchive((pPositions[i][j], pBitStrings[i][j]))

        with concurrent.futures.ThreadPoolExecutor() as executor:
            for i in range(len(f)):
                for j in range(n):
                    executor.submit(multiThreaded, i, j)
            # mapped = {executor.submit(multiThreaded, )}

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

    return sArchive


def testing(seqs, i):
    print("\nTest " + str(i))
    t = MSAMGPSO(seqs, [-2.0, 2.0], 30, 0.729844, 1.49618, 1.49618, 1.49618, .5, 3, float('inf'), 500,
                 [float('inf'), -float('inf')], 5000)
    for res in t:
        for string in bitsToStrings(res[1], seqs):
            print(string)
        print(numOfAlignedChars(bitsToStrings(res[1], seqs)))
        print(numOfInsertedIndels(res[1], seqs))


# testing(test1, 1)
# testing(test2, 2)
# testing(test3, 3)
# testing(test4, 4)
testing(test5, 5)
testing(test6, 6)
testing(test7, 7)
