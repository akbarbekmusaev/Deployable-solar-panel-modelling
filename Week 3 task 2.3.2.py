import numpy as np
from My_Euler_Solver import solve_Euler
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#import matplotlib.pyplot as plt
m = 600 #kg
P = 500 #N
c = 1 #Ns2m-2
#
def p(t):
    p = (P/2)*(1+np.cos(np.pi*t))
    return p

def f(t, v): 
    f = (p(t) - c*v**2)/m
    return f

#data = solve_Euler(f, 15, 5, 0.1)


sol = solve_ivp(f,(0, 15), [5], rtol = 0.0000001, max_step = 0.001)
tSollvp = sol.t
vSollvp = sol.y[0, :]
Euler = solve_Euler(f,15, 5, 0.000001)
vEuler = Euler[1, :] 
tEuler = Euler[0, :]
print(vEuler[-1])
print(vSollvp[-1])
plt.plot(tEuler, vEuler)
plt.plot(tSollvp,vSollvp)

maxvel = np.max(Euler, axis = 1)[1]
indices = np.where(Euler == maxvel)
print(Euler[0, indices[1]])

maxvelsvp = np.max(vSollvp)
indicesvp = np.where(vSollvp == maxvelsvp)
print(tSollvp[indicesvp])
