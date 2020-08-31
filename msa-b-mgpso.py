from typing import List
from util import numOfAlignedChars, numOfInsertedIndels
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
