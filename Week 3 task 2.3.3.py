import numpy as np
from My_Euler_Solver import solve_Euler
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#import matplotlib.pyplot as plt
m = 600 #kg
P = 500 #N
c = 1 #Ns2m-2
A = 1e-3
gr = 9.81
T = 100
v0 = 0
x0 = 0

def g(x):
    g = 2*A*x
    return g

def s(x):
    s = (g(x)*m*gr)/(np.sqrt(g(x)**2 + 1))
    return s

def p(t):
    p = P*(1-np.exp(-t/2))
    return p

def d(v):
    d = c*v*abs(v)
    return d

def f(t, z):
    difz = [z[1], 
            (p(t) - d(z[1]) - s(z[0]))/m]
    return difz


sol = solve_ivp(f, (0, T), [x0 ,v0], rtol = 0.000001)

x = sol.y[0, :]
v = sol.y[1, :]

print(sol.t[x==x.max()][0])

plt.plot(sol.t, sol.y[0, :])








