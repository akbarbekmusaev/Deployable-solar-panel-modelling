import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum
from matplotlib.animation import PillowWriter
import webbrowser

# constants
g = 9.81
l = 1
# initial conditions
x0 = 0*np.pi/6
v0 = 0
T = 4
Moment = 9.82 #Nm
m = 1 #kg


def f(t, z):
    difz = [z[1],
            (Moment/(m*l**2))-(g * np.sin(z[0])) / l]
    return difz


sol = solve_ivp(f, (0, T), [x0+np.pi/2, v0], rtol=0.00001)
x = sol.y[0, :]
v = sol.y[1, :]

plt.plot(sol.t, x)
plt.show()
ani = animate_pendulum(sol.t, x)
ani.save('animation.gif', writer=PillowWriter(fps=30))
webbrowser.open('animation.gif')