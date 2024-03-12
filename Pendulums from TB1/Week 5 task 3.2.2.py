import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_sprung_pendulum

#constants
g = 9.81
l = 0.75
x0 = np.pi/2
v0 = 0
T = 5
m = 1
Ls = 0.5
Ld = 0.6
c = 1
k = 50
rtol = 1e-6

def d(t, z):
    difz = [z[1],
            -(g/l)*np.sin(z[0])-((Ls**2*k)/(m*l**2))*z[0]-((Ld**2*c)/(m*l**2))*z[1]]
    return difz

def f(t, z):
    difz = [z[1], 
            -(g*np.sin(z[0]))/l]
    return difz

def my_eventf(t, z):
    return z[0]
my_eventf.terminal = True
my_eventf.direction = -1

def my_eventd(t, z):
    return z[0]
my_eventd.terminal = True
my_eventd.direction = +1

x = np.array([])
v = np.array([])
t = np.array([])
x = np.append(x, x0)
v = np.append(v, v0)
t = np.append(t, 0)

while (t[-1] < 20):
    sol1 = solve_ivp(f, (t[-1], t[-1] + T), [x[-1] ,v[-1]], rtol = rtol, events = my_eventf)
    x = np.append(x, sol1.y[0, :])
    v = np.append(v, sol1.y[1, :])
    t = np.append(t, sol1.t)
    if(t[-1] > 20):
        break
    sol2 = solve_ivp(d, (t[-1], t[-1] + T), [x[-1] ,v[-1]], rtol = rtol, events = my_eventd)
    x = np.append(x, sol2.y[0, :])
    v = np.append(v, sol2.y[1, :])
    t = np.append(t, sol2.t)

plt.plot(t[t<20], x[t<20])

ani = animate_sprung_pendulum(t[t<20], x[t<20])

