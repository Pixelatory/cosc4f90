"""
    Shared Functions
    Nick Aksamit 2020

    Simply contains functions that may be shared between the MSABPSO and the MSAAMPSO
"""


def aggregatedFunction(position, seq, w1, w2, infeasible):
    """A maximization aggregated fitness function that follows the following formula:

    f(x) = w1 * numOfAlignedChars(x) + w2 * (nMax - nI),\n
    where nMax is the number of total indels,\n
    and nI is the number of indels in-between characters.

    Note: if the position vector is invalid and infeasible is true, then -float('inf') is returned


    :param position: position vector
    :type position: list of (list of int)
    :param seq: sequences to be aligned
    :type seq: list of str
    :param w1: weight coefficient for number of aligned characters
    :type w1: float
    :param w2: weight coefficient for number of leading indels used
    :type w2: float
    :param infeasible: whether or not a position can be infeasable
    :type infeasible: bool
    :rtype: float
    :return: fitness value
    """

    strings = posToStrings(position, seq)

    nMax = 0  # total number of indels

    for bitlist in position:
        for bit in bitlist:
            if bit == 1:
                nMax = nMax + 1

    nI = 0  # number of indels before last char

    # The small procedure below counts for nI,
    # and also eliminates infeasible positions.
    # If the number of 0 bits is more than the
    # amount of characters for the sequence, then
    # it's invalid. (return -infinity)
    for i in range(len(seq)):
        tmp = 0
        hitLastChar = False
        for bit in position[i]:
            if bit == 0:
                tmp = tmp + 1
                if tmp == len(seq[i]):
                    hitLastChar = True
            elif bit == 1 and not hitLastChar:
                nI = nI + 1  # an indel was found before the last character in sequence

        if infeasible and tmp != len(seq[i]):
            return float('-inf')  # return a very small number, this solution is infeasible

    return (w1 * numOfAlignedChars(strings)) + (w2 * (nMax - nI))


def posToStrings(position, seq):
    """Converts a list of sequences into a list of strings with indels, according to the position vector given.

    :type position: list of (list of int)
    :type seq: list of str
    :rtype: list of str
    """
    result = []
    i = 0
    for bitlist in position:
        j = 0
        result.append("")
        for bit in bitlist:
            if bit == 0 and j < len(seq[i]):
                result[len(result) - 1] = result[len(result) - 1] + seq[i][j]
                j = j + 1
            else:
                result[len(result) - 1] = result[len(result) - 1] + "-"
        i = i + 1
    return result


def numOfAlignedChars(strings):
    """Counts the number of aligned characters in a list of strings.

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
    for i in range(len(strings[0])):
        charList.append({})
        for string in strings:
            if string[i] != "-":
                if string[i] in charList[i]:
                    charList[i][string[i]] = charList[i][string[i]] + 1
                else:
                    charList[i][string[i]] = 1

    for d in charList:
        for v in d.values():
            if v > 1:
                result = result + v
    return result
