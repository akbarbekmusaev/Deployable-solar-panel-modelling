import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum
from matplotlib.animation import PillowWriter
import webbrowser


# constants
theta_finishing = 90*(np.pi/180)
g = 9.81
# initial conditions
theta = 30
thetaspeed = 0
T = 15
Tstall = 0.19
OmegaMAx = 5000
gearRatio = 400
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
M_total = M_solar + 2*(M_driving + M_green + M_frame)
R_cm = (M_solar * R_solar + 2 * M_driving * R_driving + 2 * M_green * R_green + 2 * M_frame * R_frame) / M_total

#Moment
Moment_of_inertia = M_total * R_cm**2

def create_animation(sol, theta_finishing):
    t_animation = sol.t
    x_animation = sol.y[0, :] + np.pi/2
    ani = animate_pendulum(t_animation, x_animation)
    ani.save('animation.gif', writer=PillowWriter(fps=30))
    webbrowser.open('animation.gif')
    
def create_animation_spyder(sol, theta_finishing):
    ani = animate_pendulum(sol.t, sol.y[0, :] +np.pi/2)
    return ani

def motortorque(thetaspeed):
    torque = -(Tstall/OmegaMAx)*(thetaspeed * gearRatio) + Tstall
    return torque
def torque(thetaspeed):
    torque = (motortorque(thetaspeed) * gearRatio)
    return torque

def create_plot(x, v, theta_finishing):
    fig, axs = plt.subplots(2)
    axs[0].plot(sol.t, x, label='Angular displacement')
    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].legend()
    axs[1].plot(sol.t, v, label='Angular velocity (rad/s)')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].legend()
    plt.show()

def f(t, z):
    difz = [z[1],
            (torque(z[1])/(M_total*R_cm**2))-(g * np.cos(z[0])) / R_cm]
    return difz

def my_eventstop(t, z):
    return z[0] - theta_finishing;
my_eventstop.terminal = True
my_eventstop.direction = 1


sol = solve_ivp(f, (0, T), [theta, thetaspeed], rtol=0.00001)
x = sol.y[0, :]
v = sol.y[1, :]
create_plot(x, v, theta_finishing)
# if sol.t_events[0].size == 0:
#     print("Time to reach the finishing angle: None")
# else:
#     print("Time to reach the finishing angle: " + sol.t_events[0][0].__str__())
print("Finishing angle: " + theta_finishing.__str__())
print("Finishing angular velocity: " + sol.y[1, -1].__str__())
print("Finishing angular displacement: " + sol.y[0, -1].__str__())
print("R_cm: " + R_cm.__str__())
print("M_total: " + M_total.__str__())
print("Moment_of_inertia: " + Moment_of_inertia.__str__())
print(M_total)
#create_animation_spyder(sol, theta_finishing)
#create_animation(sol, theta_finishing)
#kghuwrg