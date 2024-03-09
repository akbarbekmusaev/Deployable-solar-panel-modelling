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
x0 = 1
v0 = 1
T = 10


def f(t, z):
    difz = [z[1],
            -(g * np.sin(z[0])) / l]
    return difz


sol = solve_ivp(f, (0, T), [x0, v0], rtol=0.00001)
x = sol.y[0, :]
v = sol.y[1, :]

plt.plot(sol.t, x)
plt.show()
ani = animate_pendulum(sol.t, x)


# Save the animation as a GIF
ani.save('animation.gif', writer=PillowWriter(fps=30))

# Open the GIF in the default web browser
webbrowser.open('animation.gif')