
import random

def optimizeBMGPSO(seq, n, w, c1, c2, bounds, maxIter):
    # Checking for trivial errors first
    if(len(bounds) != 2 or bounds[0] > bounds[1]):
        raise Exception("Bounds is inputted incorrectly")
    if(n < 1):
        raise Exception("Swarm size cannot be < 1")

    # initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []
    gBest = 0  # indice of the global best particle
    sArchive = []

    # Begin with initializing the swarm
    for i in range(n):
        position = []
        velocity = []

        # pos[0] -> n
        position.append(random.randint(0, 50))
        velocity.append(0)

        # pos[1] -> w
        wTmp = random.random() if random.uniform(0, 1) > 0.5 else -random.random()
        position.append(wTmp)
        velocity.append(0)

        cTmp = random.uniform(0, order2Stability(wTmp))

        # pos[2] -> c1
        # pos[3] -> c2
        # pos[4] -> c3
        for i in range(3):
            position.append(cTmp)
            velocity.append(0)

        # pos[5] -> maxIter
        position.append(random.randint(0, 10000))
        velocity.append(0)


        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)
        if f(position) < f(pPositions[gBest]):
            gBest = i

    # Now this is where the iterations begin
    it = 1
    while it < maxIter and f(pPositions[gBest]) > term:
        for i in range(n): # for each particle
            r1 = random.random()
            r2 = random.random()
            for j in range(dim): # update the velocity and particle position using the formulas
                pVelocities[i][j] = w * pVelocities[i][j] + \
                                    r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) + \
                                    r2 * c2 * (pPositions[gBest][j] - pPositions[i][j])
                pPositions[i][j] = pPositions[i][j] + pVelocities[i][j]

            if f(pPositions[i]) < f(pPersonalBests[i]): # update personal best if applicable
                pPersonalBests[i] = pPositions[i]
            if f(pPositions[i]) < f(pPositions[gBest]): # update global best if applicable
                gBest = i
        it = it + 1
    return pPositions[gBest]

def order2Stability(w):
    """ Returns the result of the RHS order-2 stability equation.

    :type w: float
    :param w: inertia coefficient
    :rtype: float
    """

    # |w| < 1 & 0 < c must be satisfied for stability
    if abs(w) >= 1:
        raise Exception("|w| >= 1")

    return (12 * (1 - w ** 2)) / (7 * (w + 1) - 12 * w)
def fitness(sArchive):
