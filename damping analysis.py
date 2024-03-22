from Initial_model import openingmodel
from Initial_model import closingmodel
from Initial_model import dampertorque
import matplotlib.pyplot as plt
import numpy as np

# Constants
theta_f = 90 * (np.pi / 180)
theta_i = np.pi / 12
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gear_ratio = 270
c_values = [100]
k = 50

def rad_s_to_rpm(rad_s):
    return (rad_s * 60) / (2 * np.pi)

# Define function to plot time and position of model
def TimeAndPositionOfModelClosing_plot(theta_initial, theta_finishing, c_values):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].set_title('Angular velocity of the model over Time')

    for c in c_values:
        sol = closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c)
        axs[0].plot(sol.t, sol.y[0, :], label='c = {}'.format(c))
        axs[1].plot(sol.t, sol.y[1, :], label='c = {}'.format(c))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def TimeAndPositionOfModelOpening_plot(theta_initial, theta_finishing, c_values):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].set_title('Angular velocity of the model over Time')

    for c in c_values:
        sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k)
        axs[0].plot(sol.t, sol.y[0, :], label='c = {}'.format(c))
        axs[1].plot(sol.t, sol.y[1, :], label='c = {}'.format(c))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

TimeAndPositionOfModelOpening_plot(theta_i, theta_f, c_values)
#TimeAndPositionOfModelClosing_plot(theta_f, theta_i, c_values)
# Run the opening model
sol = closingmodel(theta_f, theta_i, speed_initial, T_stall, omega_max, gear_ratio, c_values[0])

# Extract theta and speed values from the solution
theta_values = sol.y[0]
speed_values = sol.y[1]

# Calculate the damper torque for each pair of theta and speed values
torque_values = [dampertorque(c_values[0], theta, speed) for theta, speed in zip(theta_values, speed_values)]

# Plot damper torque vs theta
plt.figure()
plt.plot(theta_values, torque_values)
plt.xlabel('Theta (rad)')
plt.ylabel('Damper Torque (Nm)')
plt.title('Damper Torque vs Theta')
#plt.show()
