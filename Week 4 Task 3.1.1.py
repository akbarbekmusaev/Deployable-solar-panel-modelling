import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum


#constants
g = 9.81
l = 0.75
x0 = 1
v0 = 1
T = 10

def f(t, z):
    difz = [z[1], 
            -(g*np.sin(z[0]))/l]
    return difz


sol = solve_ivp(f, (0, T), [x0 ,v0], rtol = 0.00001)
solve_ivp

x = sol.y[0, :]
v = sol.y[1, :]

plt.plot(sol.t, x)

#ani = animate_pendulum(sol.t, x)

