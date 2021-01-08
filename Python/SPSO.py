import math
import random

'''
    The SPSO algorithm.
    
    n -> swarm size (integer)
    dim -> number of dimensions (integer)
    w -> inertia coefficient (float)
    c1 -> cognitive coefficient (float)
    c2 -> social coefficient (float)
    bounds -> upper and lower limit of search area ([Lower Limit, Upper Limit])
    term -> termination criteria (float)
    maxIter -> maximum iteration limit (integer > 0)
    f -> the function which tests for fitness (function)
    
    Initialization Process:
    Particle positions -> each xi = U(lower bound, upper bound), where U is uniformly distributed random value
    Particle personal best -> Same as initialized position
    Particle velocities -> each vi = 0
    
    PSO Topology: Star (gBest)
    Velocity Update: w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi), where r1 and r2 are uniformly random values between [0,1]
    Position Update: x(t+1) = x(t) + v(t), where t denotes the iteration number
'''
def SPSO(n,dim,w,c1,c2,bounds,term,maxIter,f):
    # Checking for trivial errors first
    if(len(bounds) != 2 or bounds[0] > bounds[1]):
        raise Exception("Bounds is inputted incorrectly")
    if(n < 1):
        raise Exception("Swarm size cannot be < 1")
    if(dim < 1):
        raise Exception("Maximum dimension number cannot be < 1")
    if not callable(f):
        raise Exception("The f parameter is not callable")

    # initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []
    gBest = 0 # indice of the global best particle

    # Begin with initializing the swarm
    for i in range(n):
        position = []
        velocity = []
        for j in range(dim):
            position.append(random.uniform(bounds[0],bounds[1]))
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


'''
    The PSOs.BPSO algorithm.

    n -> swarm size (integer > 0)
    dim -> number of dimensions (integer > 0)
    w -> inertia coefficient (float)
    c1 -> cognitive coefficient (float)
    c2 -> social coefficient (float)
    term -> termination criteria (float)
    maxIter -> maximum iteration limit (integer > 0)
    f -> the function which tests for fitness (function)

    Initialization Process:
    Particle positions -> each xi = U(0, 1), where U is uniformly distributed random value that's either 0 or 1
    Particle personal best -> Same as initialized position
    Particle velocities -> each vi = 0

    PSO Topology: Star (gBest)
    Velocity Update: w * vi + r1 * c1 * (yi - xi) + r2 * c2 * (ŷi - xi), where r1 and r2 are uniformly random values between [0,1]
    Probability on Velocity = p(t) = Sigmoid(v(t+1)), where v(t+1) is the new velocity vector
    Position Update: x(t+1) = {
                                1, if U(0,1) < p(t)
                                0, otherwise
                              }
'''
def PSOs.BPSO(n, dim, w, c1, c2, term, maxIter, f):
    # Checking for trivial errors first
    if (n < 1):
        raise Exception("Swarm size cannot be < 1")
    if (dim < 1):
        raise Exception("Maximum dimension number cannot be < 1")
    if not callable(f):
        raise Exception("The f parameter is not callable")

    # initialize the data containers
    pPositions = []
    pPersonalBests = []
    pVelocities = []
    gBest = 0  # indice of the global best particle

    # Begin with initializing the swarm
    for i in range(n):
        position = []
        velocity = []
        for j in range(dim):
            position.append(random.randint(0,1))
            velocity.append(0)
        pPositions.append(position)
        pPersonalBests.append(position)
        pVelocities.append(velocity)
        if f(position) < f(pPositions[gBest]):
            gBest = i

    # Now this is where the iterations begin
    it = 1
    while it < maxIter and f(pPositions[gBest]) > term:
        for i in range(n):  # for each particle
            r1 = random.random()
            r2 = random.random()
            for j in range(dim):  # update the velocity and particle position using the formulas
                pVelocities[i][j] = w * pVelocities[i][j] + \
                                    r1 * c1 * (pPersonalBests[i][j] - pPositions[i][j]) + \
                                    r2 * c2 * (pPositions[gBest][j] - pPositions[i][j])
                probability = Sigmoid(pVelocities[i][j])
                pPositions[i][j] = 1 if random.uniform(0,1) < probability else 0

            if f(pPositions[i]) < f(pPersonalBests[i]):  # update personal best if applicable
                pPersonalBests[i] = pPositions[i]
            if f(pPositions[i]) < f(pPositions[gBest]):  # update global best if applicable
                gBest = i
        it = it + 1
    return pPositions[gBest]

def Rastrigin(x):
    sum = 10 * len(x)
    for xval in x:
        sum += (xval ** 2) - (10 * math.cos(2 * math.pi * xval))
    return sum

def Sigmoid(x):
    return 1 / (1 + math.exp(-x))

print(Sigmoid(0))

'''
best = PSO(50,2,0.7298,1.49618,1.49618,[-5.12,5.12],0.5,5000,Rastrigin)
for i in range(15):
    result = PSO(50,2,1,0.5,0.5,[-5.12,5.12],0.5,5000,Rastrigin)
    if Rastrigin(best) > Rastrigin(result):
        best = result
    print("Done " + str(i+1) + "/15")

print(best, Rastrigin(best))'''