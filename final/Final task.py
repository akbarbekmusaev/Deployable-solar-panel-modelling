# Import necessary functions and libraries
from Trebuchet_Functions import animate_trebuchet
from Trebuchet_Functions import trebuchet_demo
import numpy as np
import matplotlib.pyplot as plt
from Finalfunc import finalfunc


# Set constants and initial parameters
Hlim = 10
stepsize1 = 0.5
stepsize2 = 0.1
res = 0.1
p = [9.81, 12, 5000, 200, 12, 9.5, 8.5, 2.0, 10]


# Function to check physical meaning of parameters
def physicalmeaningcheck(h, s, bc):
    return p[4] - bc - res > h and h <= Hlim and h > bc + 0.1 and s < p[4]  - bc + h-res

# Function to iterate through parameter space and find valid combinations
def iteratingf(H_loop, S_loop, BC_loop):
    fig, ax1 = plt.subplots()
    ax1.set_aspect('equal', adjustable='box')
    finalx = np.array([])
    params = []
    for h in H_loop:
        for bc in BC_loop:
            for s in S_loop:
                if physicalmeaningcheck(h, s, bc):
                    S1, S2, S3 = finalfunc(h, s, bc)
                    if S3.y[0, -1] > 0 and S3.y[1, -1] > 0:
                        finalx = np.append(finalx, S3.y[0, -1])
                        params.append([[h, s, bc]])
                        ax1.plot(S3.y[0, :], S3.y[1, :])
    return finalx, params

# Function to post-process and find the parameters corresponding to the maximum value
def postprocess(finalx, params):
    fig, ax2 = plt.subplots()
    ax2.plot(finalx, range(0, len(params)))
    iMax = np.argmax(finalx)
    print(max(finalx))
    print(params[iMax][0])
    return params[iMax]

# Function to create animation
def anim(p):
    S1, S2, S3 = finalfunc(p[5], p[6], p[7])
    ani = animate_trebuchet(S1.t, S1.y, S2.t, S2.y, S3.t, S3.y, p)
    return ani

# Function to find parameters for maximum distance
def findmaxdistpar(Hlim):
    H1 = np.arange(6, Hlim, stepsize1)
    S1 = np.arange(6, 10, stepsize1)
    BC1 = np.arange(1, 5, stepsize1)
    finalx1, params1 = iteratingf(H1, S1, BC1)
    finalpar1 = postprocess(finalx1, params1)[0]

    H2 = np.arange(finalpar1[0] - 0.9, finalpar1[0] + 0.9, stepsize2)
    S2 = np.arange(finalpar1[1] - 0.9, finalpar1[1] + 0.9, stepsize2)
    BC2 = np.arange(finalpar1[2] - 0.6, finalpar1[2] + 0.9, stepsize2)
    finalx2, params2 = iteratingf(H2, S2, BC2)
    finalpar2 = postprocess(finalx2, params2)[0]
    return finalpar2

# Function to find maximum error for a given maximum distance and parameters
def findmaxerror(maxdist, Hlim, p):
    S_loop = np.arange(5, 9.5, stepsize2)
    dist_loop = np.arange(50, maxdist, 1)
    error = []
    for dist in dist_loop:
        disterr = []
        for S in S_loop:
            if physicalmeaningcheck(p[5], S, p[7]):
                S1, S2, S3 = finalfunc(p[5], S, p[7])
                if S3.y[0, -1] > 50:
                    disterr.append(abs((dist - S3.y[0, -1])/dist))
        error.append(min(disterr))
    return max(error)

# Uncomment the line below if you want to find parameters for maximum distance
finalparams = findmaxdistpar(Hlim)

#maxdist = 347.6

# Uncomment the line below if you want to find maximum error for a given distance
#error = findmaxerror(maxdist, Hlim, p)

# Uncomment the line below if you want to create an animation
#ani = anim(p)
