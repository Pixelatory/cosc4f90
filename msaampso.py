import random
import math
import logging
import os
import datetime
import copy

'''
    AMPSO for the MSA Problem
    Nick Aksamit 2020

    Acknowledgement goes towards: 
'''


def MSAAMPSO(seq, genRange, coefLimits, n, w, c1, c2, vmax, vmaxiterlimit, term, maxIter, f, w1, w2, log):
    """
    Stuff here
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
    elif coefLimits != type(list) and coefLimits is not None:
        raise Exception("Coefficient limit must be a list or none")
    elif len(coefLimits) < 1:
        raise Exception("coefLimits list must be of size 1 or 2")
    elif genRange != type(list) and genRange is not None:
        raise Exception("genRange must be a list or none")
    elif len(genRange) < 1:
        raise Exception("genRange list must be of size 1 or 2")

    # sort coefLimits and genRange, just makes life easier
    if len(coefLimits) >= 2:
        if coefLimits[0] > coefLimits[1]:
            tmp = coefLimits[0]
            coefLimits[0] = coefLimits[1]
            coefLimits[1] = tmp

    if len(genRange) >= 2:
        if genRange[0] > genRange[1]:
            tmp = genRange[0]
            genRange[0] = genRange[1]
            genRange[1] = tmp

    if log:
        mkdir("MSAAMPSO logs")
        logging.basicConfig(filename=("MSAAMPSO logs/msaampso" + str(datetime.datetime.now().time()) + ".log"),
                            format='%(message)s',
                            level=logging.INFO)
        logging.info("genRange: " + str(genRange))
        logging.info("coefLimits: " + str(coefLimits))
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
        of the angular modulation formula. Then, sample random values within genRange with
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
                val = gen(random.uniform(genRange[0], genRange[1]), pos[0], pos[1], pos[2], pos[3])
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
    colLength = math.ceil(lSeq["len"] * 1.2)

    gBestPos = [0, 0, 0, 0]

    # Initializing the particles of swarm
    for i in range(n):
        position = []
        velocity = []

        # The angular modulation coefficients are the position
        for j in range(4):
            velocity.append(0)
            position.append(random.uniform(coefLimits[0], coefLimits[1]))

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


def mkdir(path):
    os.mkdir(path)
