import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum
from matplotlib.animation import PillowWriter
import webbrowser

# constants
theta_finishing = np.pi/18
g = 9.81
# initial conditions
theta_initial = np.pi/2
speed_intial = 0
T = 10
T_stall = 0.19
omega_max = 5000
gearRatio = 1000

# Length and mass constants of base components
L_base = 0.1  # m
L_longbeam = 0.8  # m
L_shortbeam = 0.4  # m
L_solarpanel = 0.67  # m
M_solarpanel = 11  # kg
rho_beams = 0.71  # kg/m

# Masses of components
M_panels = 3 * M_solarpanel
M_driving = rho_beams * L_shortbeam
M_green = rho_beams * L_shortbeam
M_blue = rho_beams * L_shortbeam
M_frame = 2 * rho_beams * L_longbeam + 4 * rho_beams * L_shortbeam
M_total = M_panels + 2 * (M_driving + M_green + M_frame + M_blue)


# Define function to calculate phi
def phi(theta):
    phi = theta + np.arcsin(L_base / L_shortbeam)
    return phi

# Define function to calculate centre of mass
def centre_of_mass(theta):
    # Calculate CoM for each component H is the horizontal distance from the pivot to the CoM
    # and V is the vertical distance

    # CoM of solar panels
    R_solar = np.sqrt((0.5 * L_solarpanel) ** 2 + L_longbeam ** 2 - L_solarpanel * L_longbeam * np.cos(phi(theta)))
    H_solar = L_longbeam * np.cos(theta) - 0.5 * L_solarpanel * np.cos(phi(theta) - theta)
    V_solar = np.sqrt(R_solar ** 2 - H_solar ** 2)
    angle_solar = np.arctan(V_solar / H_solar) * 180 / np.pi

    # CoM of top frame
    R_frame = np.sqrt((0.5 * L_longbeam) ** 2 + L_longbeam ** 2 - L_longbeam ** 2 * np.cos(phi(theta)))
    H_frame = L_longbeam * np.cos(theta) - 0.5 * L_longbeam * np.cos(phi(theta) - theta)
    V_frame = np.sqrt(R_frame ** 2 - H_frame ** 2)
    angle_frame = np.arctan(V_frame / H_frame) * 180 / np.pi

    # CoM of driving beam
    R_driving = 0.25 * L_longbeam
    H_driving = 0.25 * L_longbeam * np.cos(theta)
    V_driving = np.sqrt(R_driving ** 2 - H_driving ** 2)
    angle_driving = np.arctan(V_driving / H_driving) * 180 / np.pi

    # CoM of green beam
    R_green = np.sqrt((0.5 * L_longbeam) ** 2 + (0.25*L_longbeam) ** 2 - 0.25 * L_longbeam ** 2 * np.cos(np.pi - phi(theta)))
    H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
    V_green = np.sqrt(R_green ** 2 - H_green ** 2)
    angle_green = np.arctan(V_green / H_green) * 180 / np.pi

    #CoM of blue beam
    R_blue = np.sqrt(L_shortbeam**2 + (0.5 * L_longbeam) ** 2 - L_shortbeam ** 2 * np.cos(np.pi - phi(theta)))
    H_blue = np.sqrt(L_shortbeam ** 2 - L_base ** 2) + 0.5 * L_shortbeam * np.cos(theta)
    V_blue = np.sqrt(R_blue ** 2 - H_blue ** 2)
    angle_blue = np.arctan(V_blue / H_blue) * 180 / np.pi

    # H centre of mass
    H_cm = (M_panels * H_solar + M_frame * H_frame + M_driving * H_driving + M_green * H_green + M_blue * H_blue) / M_total

    # V centre of mass
    V_cm = (M_panels * V_solar + M_frame * V_frame + M_driving * V_driving + M_green * V_green + M_blue * V_blue) / M_total

    # R centre of mass
    R_cm = np.sqrt(H_cm ** 2 + V_cm ** 2)
    angle_cm = np.arctan(V_cm / H_cm) * 180 / np.pi

    return R_cm

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

# Define function to calculate output torque from gearbox to mechanism
def torque_out(speed):
    torque = torque_in(speed) * gearRatio
    return torque

# Define function to calculate input torque from motor to gearbox
def torque_in(speed):
    torque = -(T_stall / omega_max) * (speed * gearRatio) + T_stall
    return torque

# Define function to plot time and position of model
def TimeAndPositionOfModel_plot(x, v, theta_finishing):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].plot(sol.t, x, label='Angular displacement')
    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].legend()
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].plot(sol.t, v, label='Angular velocity (rad/s)')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].legend()
    axs[1].set_title('Angular velocity of the model over Time')
    plt.subplots_adjust(hspace=0.5)
    plt.show()

# Define function to plot change of CoM with theta
def plot_com_change():
    # Create a range of theta values from 0 to pi with a step of 0.01
    theta_values = np.arange(0, np.pi, 0.01)

    # Calculate the CoM for each theta value
    com_values = [centre_of_mass(theta) for theta in theta_values]

    # Plot the CoM values against theta values
    plt.plot(theta_values, com_values)
    plt.xlabel('Theta (rad)')
    plt.ylabel('Centre of Mass (m)')
    plt.title('Change of Centre of Mass over angle Theta')
    plt.show()

# Define function to plot motor values
def create_motor_plot(sol):
    torque_values = [torque_in(speed) for speed in sol.y[1, :]]
    motor_speed_values = [speed * gearRatio for speed in sol.y[1, :]]
    power_values = [torque * (speed * (2 * np.pi / 60)) for torque, speed in zip(torque_values, motor_speed_values)]
    total_power = np.trapz(power_values, sol.t)
    print("Total power used during the operation (J): ", total_power)
    fig, axs = plt.subplots(3, figsize=(10, 10))
    axs[0].plot(sol.t, torque_values, label='Torque (Nm)')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Torque (Nm)')
    axs[0].legend()
    axs[0].set_title('Motor torque over time')
    axs[1].plot(sol.t, motor_speed_values, label='Speed (rpm)')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Speed (rpm)')
    axs[1].legend()
    axs[1].set_title('Motor speed over time')
    axs[2].plot(sol.t, power_values, label='Power (W)')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_ylabel('Power')
    axs[2].legend()
    axs[2].set_title('Motor power over time')
    plt.subplots_adjust(hspace = 1)
    plt.show()

# Define function for differential equation
def f(t, z):
    difz = [z[1],
            (-torque_out(z[1]) / (M_total * centre_of_mass(z[0]) ** 2)) - (g * np.cos(z[0])) / centre_of_mass(z[0])]
    return difz

# Define event to stop the simulation when the pendulum reaches the finishing angle
def my_eventstop(t, z):
    return z[0] - theta_finishing;

my_eventstop.terminal = True
my_eventstop.direction = -1

# Solve the differential equation
sol = solve_ivp(f, (0, T), [theta_initial, speed_intial], rtol=0.000001, events=(my_eventstop))

#Extract the results
x = sol.y[0, :]
v = sol.y[1, :]

# Plot the position and speed of the model against time
TimeAndPositionOfModel_plot(x, v, theta_finishing)

# Plot the motor values
create_motor_plot(sol)

# Plot the change of CoM with theta
#plot_com_change()

# Print the results
if sol.t_events[0].size == 0:
    print("Time to reach the finishing angle (s): Did not reach the angle")
else:
    print("Time to reach the finishing angle (s): " + sol.t_events[0][0].__str__())
print("Finishing angular velocity of the model (rad/s): " + sol.y[1, -1].__str__())
print("Finishing angular displacement of the model (rad): " + sol.y[0, -1].__str__())
print("Finishing torque of motor (Nm): " + torque_in(sol.y[1, -1]).__str__())
print("Finishing angular velocity of motor (RPM): " + (sol.y[1, -1] * gearRatio).__str__())
print("Holding torque for the model (Nm): " + (M_total * g * centre_of_mass(theta_initial) * np.cos(theta_initial)).__str__())

# Uncomment to create animations
#create_animation_spyder(sol, theta_finishing)
create_animation(sol, theta_finishing)
