# Import necessary functions and libraries
from Trebuchet_Functions import trebuchet_demo
from Trebuchet_Functions import animate_trebuchet
import numpy as np
import matplotlib.pyplot as plt

# Function to check physical meaning of parameters
def physicalmeaningcheck(h, s, bc):
    return B - bc - res > h and h <= 10 and h > bc + 0.1 and s < B - bc + h - res

# Function to iterate through parameter space and find valid combinations
def iteratingf(H_loop, S_loop, BC_loop):
    finalx = np.array([])
    params = []
    for h in H_loop:
        for bc in BC_loop:
            for s in S_loop:
                if physicalmeaningcheck(h, s, bc):
                    S1, S2, S3 = trebuchet_demo(h, s, bc)
                    if S3.y[0, -1] > 0 and S3.y[1, -1] > 0:
                        finalx = np.append(finalx, S3.y[0, -1])
                        params.append([[h, s, bc]])
    return finalx, params

# Function to post-process and find the parameters corresponding to the maximum value
def postprocess(finalx, params):
    iMax = np.argmax(finalx)
    print(max(finalx))
    print(params[iMax][0])
    return params[iMax]

# Function to create animation
def anim(p):
    S1, S2, S3 = trebuchet_demo(p[5], p[6], p[7])
    ani = animate_trebuchet(S1.t, S1.y, S2.t, S2.y, S3.t, S3.y, p)
    return ani

# Function to find parameters for maximum distance
def findmaxdistpar(Hlim, B):
    H1 = np.arange(res, Hlim, stepsize1)
    BC1 = np.arange(res, B, stepsize1)
    S1 = np.arange(res, Hlim + B, stepsize1)
    finalx1, params1 = iteratingf(H1, S1, BC1)
    finalpar1 = postprocess(finalx1, params1)[0]

    H2 = np.arange(finalpar1[0] - 1, finalpar1[0] + 1, stepsize2)
    S2 = np.arange(finalpar1[1] - 1, finalpar1[1] + 1, stepsize2)
    BC2 = np.arange(finalpar1[2] - 1, finalpar1[2] + 1 + res, stepsize2)
    finalx2, params2 = iteratingf(H2, S2, BC2)
    finalpar2 = postprocess(finalx2, params2)[0]
    return finalpar2

# Function to find string length for a given distance and parameters
def findstringlength(dist, Hlim, p):
    S_loop = np.arange(res, Hlim + p[4], stepsize2)
    diff = []
    for S in S_loop:
        if physicalmeaningcheck(p[5], S, p[7]):
            S1, S2, S3 = trebuchet_demo(p[5], S, p[7])
            diff.append(abs(dist - S3.y[0, -1]))
    return S_loop[diff.index(min(diff))], (min(diff)/dist)*100

# Function to find maximum error for a given maximum distance and parameters
def findmaxerror(maxdist, Hlim, p):
    S_loop = np.arange(res, Hlim + p[4], stepsize2)
    dist_loop = np.arange(50, maxdist, 1)
    error = []
    for S in S_loop:
        for dist in dist_loop:
            if physicalmeaningcheck(p[5], S, p[7]):
                S1, S2, S3 = trebuchet_demo(p[5], S, p[7])
                error.append(abs(dist - S3.y[0, -1])/dist)
    return max(error)*100

# Set constants and initial parameters
B = 8
Hlim = 10
stepsize1 = 0.7
stepsize2 = 0.1
res = 0.1

# Uncomment the line below if you want to find parameters for maximum distance
# finalparams = findmaxdistpar(Hlim, B)

# Define initial parameters
H = 6.3
S = 6.1
BC = 1.6
BP = B - BC

p = [9.81, 12, 4000, 100, B, H, S, BC, BP]
maxdist = 60

# Uncomment the line below if you want to find maximum error for a given distance
# error = findmaxerror(maxdist, Hlim, p)

# Uncomment the line below if you want to create an animation
# ani = anim(p)
