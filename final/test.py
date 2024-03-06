# Import necessary functions and libraries
from Trebuchet_Functions import animate_trebuchet
from Trebuchet_Functions import trebuchet_demo
import numpy as np
import matplotlib.pyplot as plt
from Finalfunc import finalfunc


# Set constants and initial parameters
Hlim = 10
stepsize1 = 0.7
stepsize2 = 0.1
res = 0.1
p = [9.81, 12, 5000, 200, 12, 9.9, 9.5, 2.0, 10]

# Function to create animation
S1, S2, S3 = finalfunc(p[5], p[6], p[7])
ani = animate_trebuchet(S1.t, S1.y, S2.t, S2.y, S3.t, S3.y, p)
