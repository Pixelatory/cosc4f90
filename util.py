import copy
import random
import math
import os
from typing import Dict, List, Tuple

"""
    Utility
    Nick Aksamit 2020

    Simply contains functions that may be shared between the code pieces.
"""


def aggregatedFunction(bitmatrix, seq, w1, w2, checkInfeasibility=False, op=None):
    """A maximization aggregated fitness function that follows the following formula:

    f(x) = w1 * numOfAlignedChars(x) + w2 * (nMax - nI),\n
    where nMax is the number of total indels,\n
    and nI is the number of indels in-between characters.

    Note: if the position vector is invalid and infeasible is true, then -float('inf') is returned


    :param bitmatrix: position vector
    :type bitmatrix: List[List[int]]
    :param seq: sequences to be aligned
    :type seq: List[str]
    :param w1: weight coefficient for number of aligned characters
    :type w1: float
    :param w2: weight coefficient for number of leading indels used
    :type w2: float
    :param checkInfeasibility: whether or not a position can be infeasable
    :type checkInfeasibility: bool
    :param op: how to compare between count of 0 bits in row, and the length of the sequence
    :type op: (int, int) -> bool
    :rtype: float
    :return: fitness value
    """

    if checkInfeasibility and op is None:
        raise Exception("Cannot check infeasibility with a None operator")

    if checkInfeasibility and infeasible(bitmatrix, seq, op):
        return -float('inf')

    strings = bitsToStrings(bitmatrix, seq)

    nMax = maxNumOfIndels(bitmatrix, seq)  # total number of indels

    nI = numOfInsertedIndels(bitmatrix, seq)  # number of indels before last char

    return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI))


def numOfInsertedIndels(bitmatrix, seq):
    """Counts the number of indels that are found before the last character in a sequence.

    :type bitmatrix: List[List[int]]
    :type seq: List[str]
    :rtype: int
    """
    # Remember: bit 0 means character from sequence
    #           bit 1 means inserting indel
    count = 0

    for i in range(len(seq)):  # loop through each sequence
        tmp = 0
        hitLastChar = False  # whether the last char in the sequence was hit
        for bit in bitmatrix[i]:
            if bit == 0:
                tmp += 1
                if tmp == len(seq[i]):
                    hitLastChar = True
            else:
                count += 1

            if hitLastChar:
                break

    return count


def numOfIndels(bitmatrix):
    """Counts the total number of indels according to the individual bits.

    :type bitmatrix: List[List[int]]
    :rtype: int
    """
    count = 0
    for bitlist in bitmatrix:
        for bit in bitlist:
            if bit == 1:
                count += 1
    return count


def maxNumOfIndels(bitmatrix, seq):
    """Counts the total number of indels according to the difference between
    bitmatrix row length and sequence length.

    Does not consider individual bits.

    :type bitmatrix: List[List[int]]
    :type seq: List[str]
    :rtype: int
    """
    count = 0
    for i in range(len(bitmatrix)):
        count += len(bitmatrix[i]) - len(seq[i])
    return count


def infeasible(bitmatrix, seq, op):
    """Determines if a bitmatrix and sequence combination is infeasible.

    Counts the number of 0 bits in each row of the bitmatrix, and if it
    doesn't match up with an operator comparison between the 0 bit count
    and the length of the sequence, returns True for infeasibility.

    The operator can be anything within the parameter format, however
    usually it is >, <, =. For easy access to these use the operator
    module.

    Ex. When using >, result is True when [0's count] > [length of string]

    :type op: (int, int) -> bool
    :param op: comparison of (count of 0 bits in row, length of sequence)
    :type bitmatrix: List[List[int]]
    :type seq: List[str]

    """
    for i in range(len(seq)):
        count = 0
        for bit in bitmatrix[i]:
            if bit == 0:
                count += 1

        if op(count, len(seq[i])):
            return True


def bitsToStrings(bitmatrix, seq):
    """Converts a list of sequences into a list of strings with indels, according to the bitmatrix provided.

    Note: A bit of 0 means a character, and a bit of 1 means indel.
    However, if the number of 0 bits exceeds the length of the sequence,
    indels will be inserted instead.

    :type bitmatrix: List[List[int]]
    :type seq: List[str]
    :rtype: List[str]
    """
    result = []
    i = 0
    for bitlist in bitmatrix:
        j = 0
        result.append("")
        for bit in bitlist:
            if bit == 0 and j < len(seq[i]):
                result[len(result) - 1] += seq[i][j]
                j += 1
            else:
                result[len(result) - 1] += "-"
        i = i + 1
    return result


def numOfAlignedChars(strings):
    """Counts the number of aligned characters in a list of strings.

    **Assumes that each string is of the same length.**

    :type strings: List[str]
    :rtype: int
    """
    if len(strings) < 1:
        raise Exception("There's no strings in the num of aligned chars function")
    elif len(strings) == 1:
        print("Warning: only 1 string in the numOfAlignedChars function")
        return 0

    result = 0
    charList = []

    # The charList has a dictionary, and each dictionary in the list represents
    # a column in each string
    for i in range(len(strings[0])):
        charList.append({})
        for string in strings:
            if string[i] != "-":
                if string[i] in charList[i]:
                    charList[i][string[i]] += 1
                else:
                    charList[i][string[i]] = 1

    for d in charList:
        for v in d.values():
            if v > 1:
                result = result + v
    return result


def getLongestSeqDict(seq):
    """ Returns a dictionary of information on the longest sequence within a list of sequences.

    -> "idx": index of the longest sequence

    -> "len": length of longest sequence

    :rtype: Dict[str,int]
    :type seq: List[str]
    """
    lSeq = {  # longest sequence (index value and length)
        "idx": 0,
        "len": len(seq[0])
    }

    for i in range(len(seq)):
        s = len(seq[i])

        if lSeq["len"] < s:
            lSeq["idx"] = i
            lSeq["len"] = s

    return lSeq


def genBitMatrix(pos, seq, colLength, genInterval):
    """Generates a bit matrix given the position vector and the angular modulation formula (gen function).

    :type pos: List[float]
    :type genInterval: List[float]
    :type seq: List[str]
    :type colLength: int
    :rtype: List[List[int]]
    """
    bitmatrix = []
    for li in range(len(seq)):
        bitmatrix.append([])
        for ti in range(colLength):
            val = gen(random.uniform(genInterval[0], genInterval[1]), pos[0], pos[1], pos[2], pos[3])
            if val > 0:
                bitmatrix[li].append(1)
            else:
                bitmatrix[li].append(0)
    return bitmatrix


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


def dominates(seq, bm1, bm2):
    """Checking that bit matrix 1 dominates bit matrix 2.

    Used in both MGPSO py code files.

    :type bm1: List[List[int]]
    :type bm2: List[List[int]]
    :type seq: List[str]
    """
    better = False
    # numOfAlignedChars
    strings1 = bitsToStrings(bm1, seq)
    strings2 = bitsToStrings(bm2, seq)
    res1 = numOfAlignedChars(strings1)
    res2 = numOfAlignedChars(strings2)
    if res1 < res2:
        return False
    elif res1 > res2:
        better = True

    # numOfInsertedIndels
    res1 = numOfInsertedIndels(bm1, seq)
    res2 = numOfInsertedIndels(bm2, seq)
    if res1 > res2:
        return False
    elif res1 < res2:
        better = True

    return better


def archiveGuide(seq, sArchive, bmidx, distidx, k):
    """Uses tournament selection where k is the number of particles to choose.

    Out of the k possible particles randomly selected, the least crowded particle wins the tournament.

    :type seq: List[str]
    :type sArchive: List[List[List[float], List[List[int]], float]] | List[List[List[List[int]], float]]
    :type bmidx: int
    :type distidx: int
    :type k: int
    :rtype: List[float]
    :returns: Position vector
    """
    updateCrowdingDistances(seq, sArchive, bmidx, distidx)

    idx = random.randint(0, len(sArchive) - 1)

    for i in range(k - 1):
        tmp = random.randint(0, len(sArchive) - 1)
        if sArchive[tmp][2] > sArchive[idx][2]:
            idx = tmp

    return sArchive[idx][0]


def updateCrowdingDistances(seq, sArchive, bmidx, distidx):
    """Calculates the crowding distance of each archive solution.

    :type seq: List[str]
    :type bmidx: int
    :type distidx: int
    :type sArchive: List[List[List[float], List[List[int]], float]] | List[List[List[List[int]], float]]
    """
    for i in range(2):
        if i == 0:  # numOfAlignedChars
            sArchive.sort(key=lambda x: -numOfAlignedChars(bitsToStrings(x[bmidx], seq)))
        else:  # numOfInsertedIndels
            sArchive.sort(key=lambda x: numOfInsertedIndels(x[bmidx], seq))

        # set first and last to infinite distance
        sArchive[0][distidx] = sArchive[len(sArchive) - 1][distidx] = float('inf')

        for j in range(1, len(sArchive) - 2):
            if sArchive[j][distidx] != float('inf'):
                sArchive[j][distidx] += sArchive[j + 1][distidx] - sArchive[j - 1][distidx]


def theSame(x, y):
    """Checks if two sets of bit matrix values are the same.

    Assumes that the bit matrixes are of the same length.

    :type x: List[List[int]]
    :type y: List[List[int]]
    :rtype: bool
    """
    for i in range(len(x)):
        for j in range(len(x[i])):
            if x[i][j] != y[i][j]:
                return False

    return True


def addToArchive(seq, sArchive, x, bmidx, distidx, archiveLimit):
    """Adds solution x to archive a if x is not dominated by any archive solutions.

    After adding, if any particles in the archive are now dominated, then they are removed.

    Additionally, if the archive is full then the most crowded solution is removed.

    :type seq: List[str]
    :type x: Tuple[List[float], List[List[int]]] | List[List[int]]
    :param x: Either a tuple of position vector and bit matrix, or just a bitmatrix
    :type sArchive: List[List[List[float], List[List[int]], float]] | List[List[List[List[int]], float]]
    :param sArchive: A union of either one type-set of values or the other; not mixed
    :type bmidx: int
    :type distidx: int
    :type archiveLimit: int
    :param archiveLimit: Max number of particles allowed in the archive
    :rtype: List[List[List[float], List[List[int]], float]] | List[List[List[List[int]], float]]
    """
    if type(x) is list and len(seq) == len(x):  # means it's List[List[int]]
        bm = x
    elif type(x) is tuple and len(x) == 2:  # means it's Tuple[List[float], List[List[int]]]
        bm = x[1]
    else:
        raise Exception("Invalid parameter x, must be ([float], [[int]]) or [[int]]")

    aDominated = []  # all the archive components that are dominated by x

    # First, check that this new solution dominates every archive solution,
    # and that the values (bitstring and position) aren't repeated
    for s in sArchive:
        if dominates(seq, s[bmidx], bm):
            return sArchive
        elif dominates(seq, bm, s[bmidx]):
            aDominated.append(s)

        if theSame(s[bmidx], bm):
            return sArchive

    # When the function reaches here, it's safe to say the solution dominates.
    # So, let's add it to archive.
    if type(x) is list and len(seq) == len(x):  # means it's List[List[int]]
        sArchive.append(copy.deepcopy([bm, 0.0]))
    else:  # means it's Tuple[List[float], List[List[int]]]
        sArchive.append(copy.deepcopy([x[0], x[1], 0.0]))

    # Next, after adding it remove all the archive elements that x dominates
    newArchive = [x for x in sArchive if x not in aDominated]

    if len(sArchive) > archiveLimit:
        updateCrowdingDistances(seq, newArchive, bmidx, distidx)
        removeCrowdedSolution(newArchive, distidx)

    return newArchive


def removeCrowdedSolution(sArchive, idx):
    """Removes the most crowded solution in the archive.

    Note: does not update crowding distance before removal

    :type sArchive: List[List[Any]]
    :type idx: int
    :param idx: the index where crowding distance measure is stored in archive
    """
    minIdx = 0

    for i in range(len(sArchive)):
        if sArchive[i][idx] < sArchive[minIdx][idx]:
            minIdx = i

    del sArchive[minIdx]


def mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


test1 = ["AGQYHECK", "AFGPWERKYV", "ASWIELKV"]
test2 = ["GAAAGTG", "CGACACTAGA", "CGCAGT"]
test3 = ["TCATGT", "GCGAT", "CGTTGT", "TCGATT", "AGCACTAG", "GAGTAGAC"]
test4 = ["DMHCMHDHMMDDMPM", "MMDCCDCCPCPCHPDPC"]
test5 = ["SCWIISRSWIWCICCRI", "WCSIWSWIWWISRICWI", "WSWWIWRCCISWCISI", "RRCCWSIRRCSRWS", "SWCRWSWSWIIRISWI"]
test6 = ["ATAHVVTAFIIWGSSGWWQFGIGVI", "IVISFVWQTIIIAGIIQFSHGAST"]
test7 = ["CDGAGIATDAWNFWAVDECVIYQIYI", "AEYGKYITDWCQLNWNCWKFTIDQGL", "GLFKLNYGDWYDVICINIQW",
         "FNADCDVYGENKETGLCAEFAENQWC", "IGGQQNLTFDLLCTIECWQYGI", "LEKQNCQNKNTTKFIIFLDDLV",
         "QIQGLYFLANGKAVVCKNKYTTN", "QFGAGFDKAEIENCQDTYCLFQGWEQK", "GFDWETLWWLIKFYEFTGTICCWNN",
         "GEDYWAGGVKIVGGICADKAEWKA"]
