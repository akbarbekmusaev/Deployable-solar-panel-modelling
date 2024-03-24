from Initial_model import openingmodel
from Initial_model import closingmodel
from Initial_model import dampertorque
import matplotlib.pyplot as plt
import numpy as np

# Constants
theta_f = 100 * (np.pi / 180)
theta_i = 15 * (np.pi / 180)
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gear_ratio = 270
c_values = [10]
k_values = [0]

def rad_s_to_rpm(rad_s):
    return (rad_s * 60) / (2 * np.pi)

# Define function to plot time and position of model
def TimeAndPositionOfModelClosing_plot_c(theta_initial, theta_finishing,k, c_values):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].set_title('Angular velocity of the model over Time')

    for c in c_values:
        sol = closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k)
        axs[0].plot(sol.t, sol.y[0, :], label='c = {}'.format(c))
        axs[1].plot(sol.t, sol.y[1, :], label='c = {}'.format(c))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def TimeAndPositionOfModelOpening_plot_c(theta_initial, theta_finishing, k, c_values):
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
def TimeAndPositionOfModelOpening_plot_k(theta_initial, theta_finishing, k_values, c):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].set_title('Angular velocity of the model over Time')

    for k in k_values:
        sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k)
        axs[0].plot(sol.t, sol.y[0, :], label='k = {}'.format(k))
        axs[1].plot(sol.t, sol.y[1, :], label='k = {}'.format(k))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def TimeAndPositionOfModelClosing_plot_k(theta_initial, theta_finishing, k_values, c):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_xlabel('Time')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].set_title('Angular velocity of the model over Time')

    for k in k_values:
        sol = closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k)
        axs[0].plot(sol.t, sol.y[0, :], label='k = {}'.format(k))
        axs[1].plot(sol.t, sol.y[1, :], label='k = {}'.format(k))
        print(f"Opening model stopped at time: {sol.t_events[0]}")

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

#TimeAndPositionOfModelClosing_plot_k(theta_f, theta_i, k_values, c_values[0])
#TimeAndPositionOfModelClosing_plot_c(theta_f, theta_i, k_values[0], c_values)
TimeAndPositionOfModelOpening_plot_k(theta_i, theta_f, k_values, c_values[0])
#TimeAndPositionOfModelOpening_plot_c(theta_i, theta_f, k_values[0], c_values)

