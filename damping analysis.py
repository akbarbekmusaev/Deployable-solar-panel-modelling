from Initial_model import openingmodel
from Initial_model import closingmodel
import matplotlib.pyplot as plt
import numpy as np

# Constants
theta_f = 100 * (np.pi / 180)
theta_i = 5 * (np.pi / 180)
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gear_ratio = 370
c_values_opening = [900]
c_values_closing = [2500]
efficiency = 0.73

#08-020
k_start = 0
k_finish = 0

# Define function to plot time and position of model

def TimeAndPositionOfModelOpening_plot_c(theta_initial, theta_finishing, c_values, k_start, k_finish):
    fig, axs = plt.subplots(4, figsize=(10, 20))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Angular velocity (rad/s)')
    axs[1].set_title('Angular velocity of the model over Time')
    axs[2].set_xlabel('Angle (rad)')  # New subplot for Torque vs Angle
    axs[2].set_ylabel('Power (W)')
    axs[2].set_title('Power dissipated by the damper vs time')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Torque (Nm)')
    axs[3].set_title('Torque of the damper vs time')

    for c in c_values:
        sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k_start, k_finish, efficiency)
        axs[0].plot(sol.t, sol.y[0, :], label='c = {}'.format(c))
        axs[1].plot(sol.t, sol.y[1, :], label='c = {}'.format(c))
        # Calculate the torque at each time step
        power = sol.y[1, :] * c * sol.y[1, :]
        torque = sol.y[1, :] * c
        # Plot the torque against the angle
        axs[2].plot(sol.y[0, :], power, label='c = {}'.format(c))
        axs[3].plot(sol.t, torque, label='c = {}'.format(c))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def TimeAndPositionOfModelClosing_plot_c(theta_initial, theta_finishing, c_values, k_start, k_finish):
    fig, axs = plt.subplots(4, figsize=(10, 20))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time(s)')
    axs[1].set_ylabel('Angular velocity (rad/s)')
    axs[1].set_title('Angular velocity of the model over Time')
    axs[2].set_xlabel('Angle (rad)')  # New subplot for Torque vs Angle
    axs[2].set_ylabel('Power (W)')
    axs[2].set_title('Power dissipated by the damper vs time')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Torque (Nm)')
    axs[3].set_title('Torque of the damper vs time')

    for c in c_values:
        sol = closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k_start, k_finish, efficiency)
        axs[0].plot(sol.t, sol.y[0, :], label='c = {}'.format(c))
        axs[1].plot(sol.t, sol.y[1, :], label='c = {}'.format(c))
        # Calculate the torque at each time step
        power = sol.y[1, :] * c * sol.y[1, :]
        torque = sol.y[1, :] * c
        # Plot the torque against the angle
        axs[2].plot(sol.y[0, :], power, label='c = {}'.format(c))
        axs[3].plot(sol.t, torque, label='c = {}'.format(c))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

#TimeAndPositionOfModelOpening_plot_c(theta_i, theta_f, c_values_opening, k_start, k_finish)
TimeAndPositionOfModelClosing_plot_c(theta_f, theta_i, c_values_closing, k_start, k_finish)


