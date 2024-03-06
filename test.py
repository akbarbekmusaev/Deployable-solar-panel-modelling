from joblib import Parallel, delayed
import numpy as np
from Trebuchet_Functions import trebuchet_demo

B = 8
Hlim = 10
limF = Hlim + B
stepsize = 0.1

# Calculate the maximum number of iterations based on the step size
max_iterations = int(((Hlim - stepsize) / stepsize) ** 3)

# Preallocate space for finalx
finalx = np.empty(max_iterations)

index = 0

H_loop = np.arange(0.1, Hlim, stepsize)
BC_loop = np.arange(0.1, Hlim + stepsize, stepsize)
S_loop = np.arange(0.1, limF, stepsize)

# Define the function to be parallelized
def trebuchet_simulation(h, bc, s):
    S1, S2, S3 = trebuchet_demo(h, s, bc)
    if S3.y[0, -1] > 0:
        return S3.y[0, -1]
    else:
        return None

# Use Parallel and delayed to parallelize the loop
results = Parallel(n_jobs=2)(
    delayed(trebuchet_simulation)(h, bc, s) 
    for h in H_loop
    for bc in BC_loop
    for s in S_loop
    if B - bc > h and h <= 10 and h > bc + 0.1 and s < B - bc + h
)

# Filter out None values and store the results in finalx
finalx = np.array([result for result in results if result is not None])

# You can print the length of finalx to verify the number of valid results
print(len(finalx))
