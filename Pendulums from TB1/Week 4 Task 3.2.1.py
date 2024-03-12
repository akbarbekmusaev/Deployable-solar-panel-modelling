import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum

#constants
g = 9.81
l = 0.75
x0 = np.pi-0.25
v0 = 0
T = 5

def f(t, z):
    difz = [z[1], 
            -(g*np.sin(z[0]))/l]
    return difz

def my_event(t, z) :
    return z[0]-2
my_event.terminal = True
my_event.direction = 1

sol = solve_ivp(f, (0, T), [x0 ,v0], rtol = 0.00001, events = my_event)
solve_ivp

x = sol.y[0, :]
v = sol.y[1, :]

ani = animate_pendulum(sol.t, x)

