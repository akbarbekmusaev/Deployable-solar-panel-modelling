import numpy as np
from My_Euler_Solver import solve_Euler

#import matplotlib.pyplot as plt
m = 600 #kg
P = 500 #N
c = 1 #Ns2m-2

def p(t):
    p = (P/2)*(1+np.cos(np.pi*t))
    return p

def f(t, v): 
    f = (p(t) - c*v**2)/m
    return f

data = solve_Euler(f, 15, 5, 0.1)
