import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum
from matplotlib.animation import PillowWriter
import webbrowser


# constants
theta_finishing = 120*(np.pi/180)
g = 9.81
# initial conditions
theta = 0
thetaspeed = 0
T = 6
def phi(theta):
    phi = theta - np.arcsin(100 / 500)
    return phi


# Length and mass constants of beams
L_longbeam = 1  # m
L_shortbeam = 0.5  # m
rho_beams = 0.71  # kg/m

# Masses and CoM of components
M_solar = 33
M_driving = rho_beams * L_shortbeam
M_green = rho_beams * L_shortbeam
M_blue = rho_beams * L_shortbeam
M_frame = 4 * rho_beams * L_longbeam
R_solar = 1.1 - (0.67 / 2) * np.cos(phi(theta))
R_driving = 1 / 4
R_green = 0.5 + 0.25 * np.cos(phi(theta))
R_frame = np.arctan((0.5 * np.sin(phi(theta))) / (1 - 0.5 * np.cos(phi(theta))))
R_blue = np.arctan((0.1 + 0.25 * np.cos(theta)) / (200 * np.sqrt(6) + 0.25 * np.sin(theta)))

# Calculate the total mass and CoM
M_total = M_solar + M_driving + M_green + M_frame
R_cm = (M_solar * R_solar + M_driving * R_driving + M_green * R_green + M_frame * R_frame) / M_total

#Moment
Moment = M_total * g * R_cm + 0.2

def create_animation(sol, theta_finishing):
    T_event = sol.t_events[0][0]
    t_animation = sol.t[sol.t <= T_event]
    x_animation = sol.y[0, :][sol.t <= T_event] + np.pi/2
    ani = animate_pendulum(t_animation, x_animation)
    ani.save('animation.gif', writer=PillowWriter(fps=30))
    webbrowser.open('animation.gif')
    
def create_animation_spyder(sol, theta_finishing):
    ani = animate_pendulum(sol.t, sol.y[0, :] +np.pi/2)
    return ani
    

def create_plot(x, v, theta_finishing):
    fig, axs = plt.subplots(2)
    axs[0].plot(sol.t, x, label='Theta')
    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Theta')
    axs[0].legend()
    axs[1].plot(sol.t, v, label='Velocity')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Velocity')
    axs[1].legend()
    plt.show()

def f(t, z):
    difz = [z[1],
            (Moment/(M_total*R_cm**2))-(g * np.cos(z[0])) / R_cm]
    return difz

def my_eventstop(t, z):
    return z[0] - theta_finishing;
my_eventstop.terminal = True
my_eventstop.direction = 1

sol = solve_ivp(f, (0, T), [theta, thetaspeed], rtol=0.00001, events=my_eventstop)
x = sol.y[0, :]
v = sol.y[1, :]
create_plot(x, v, theta_finishing)
create_animation_spyder(sol, theta_finishing)
create_animation(sol, theta_finishing)

