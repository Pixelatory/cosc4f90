"""
    Shared Functions
    Nick Aksamit 2020

    Simply contains functions that may be shared between the MSABPSO and the MSAAMPSO
"""


def aggregatedFunction(bitmatrix, seq, w1, w2, checkInfeasability):
    """A maximization aggregated fitness function that follows the following formula:

    f(x) = w1 * numOfAlignedChars(x) + w2 * (nMax - nI),\n
    where nMax is the number of total indels,\n
    and nI is the number of indels in-between characters.

    Note: if the position vector is invalid and infeasible is true, then -float('inf') is returned


    :param bitmatrix: position vector
    :type bitmatrix: list of (list of int)
    :param seq: sequences to be aligned
    :type seq: list of str
    :param w1: weight coefficient for number of aligned characters
    :type w1: float
    :param w2: weight coefficient for number of leading indels used
    :type w2: float
    :param checkInfeasability: whether or not a position can be infeasable
    :type checkInfeasability: bool
    :rtype: float
    :return: fitness value
    """

    if checkInfeasability and infeasible(bitmatrix, seq):
        return -float('inf')

    strings = posToStrings(bitmatrix, seq)

    nMax = maxNumOfIndels(bitmatrix, seq)  # total number of indels

    nI = numOfInsertedIndels(bitmatrix, seq)  # number of indels before last char

    return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI))


def numOfInsertedIndels(bitmatrix, seq):
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
    count = 0
    for bitlist in bitmatrix:
        for bit in bitlist:
            if bit == 1:
                count += 1
    return count


def maxNumOfIndels(bitmatrix, seq):
    count = 0
    for i in range(len(bitmatrix)):
        count += len(bitmatrix[i]) - len(seq[i])
    return count


def infeasible(bitmatrix, seq):
    for i in range(len(seq)):
        count = 0
        for bit in bitmatrix[i]:
            if bit == 0:
                count += 1

        if count < len(seq[i]):
            return True


def posToStrings(bitmatrix, seq):
    """Converts a list of sequences into a list of strings with indels, according to the bitmatrix given.

    :type bitmatrix: list of (list of int)
    :type seq: list of str
    :rtype: list of str
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

    Assumes that each string is of the same length.

    :type strings: list of str
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
