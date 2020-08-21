import random
import math
import os
from typing import Dict, List, Tuple

"""
    Utility
    Nick Aksamit 2020

    Simply contains functions that may be shared between the code pieces.
"""


def aggregatedFunction(bitmatrix, seq, w1, w2, checkInfeasability=False, op=None):
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
    :param checkInfeasability: whether or not a position can be infeasable
    :type checkInfeasability: bool
    :param op: comparison of (count of 0 bits in row, length of sequence)
    :type op: (int, int) -> bool
    :rtype: float
    :return: fitness value
    """

    if checkInfeasability and op is None:
        raise Exception("Cannot check infeasability with a None operator")

    if checkInfeasability and infeasible(bitmatrix, seq, op):
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
