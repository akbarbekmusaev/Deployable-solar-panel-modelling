from Initial_model import openingmodel
from Initial_model import closingmodel
from Pendulum_Animations import animate_pendulum
from Pendulum_Animations import animate_sprung_pendulum
from matplotlib.animation import PillowWriter
import webbrowser
import numpy as np

#constants
theta_finishing = np.pi / 18
theta_initial = np.pi / 2
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gearRatio = 270
c = 400

# Define function to create animation
def create_animation(sol, theta_finishing):
    t_animation = sol.t
    x_animation = sol.y[0, :] + np.pi / 2
    ani = animate_pendulum(t_animation, x_animation)
    ani.save('animation.gif', writer=PillowWriter(fps=30))
    webbrowser.open('animation.gif')

# Define function to create animation for Spyder
def create_animation_spyder(sol, theta_finishing):
    ani = animate_pendulum(sol.t, sol.y[0, :] + np.pi / 2)
    return ani

# Define function to create animation for Spyder
def create_animation_spyder_damper(sol, theta_finishing):
    ani = animate_sprung_pendulum(sol.t, sol.y[0, :] + np.pi / 2)
    return ani
sol = closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio, c)
create_animation(sol, theta_finishing)