from Initial_model import openingmodel
from Initial_model import closingmodel
import matplotlib.pyplot as plt

# Constants
theta_f = 100 * (np.pi / 180)
theta_i = 5 * (np.pi / 180)
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gear_ratio = 270
c_values_opening = [500]
c_values_closing = [2500]
k_start = 0
k_finish = 0

def TimeAndPositionOfModelOpening_plot_c(theta_initial, theta_finishing, c_values, k_start, k_finish):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Moment (Nm)')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_title('Moment of the damper over Time for Opening Model')

    for c in c_values:
        sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k_start, k_finish)
        # Calculate the moment at each time step
        moment = sol.y[1, :] * c
        # Plot the moment against time
        axs[0].plot(sol.t, moment, label='c = {}'.format(c))

    axs[0].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def TimeAndPositionOfModelClosing_plot_c(theta_initial, theta_finishing, c_values, k_start, k_finish):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_ylabel('Moment (Nm)')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_title('Moment of the damper over Time for Closing Model')

    for c in c_values:
        sol = closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, k_start, k_finish)
        # Calculate the moment at each time step
        moment = sol.y[1, :] * c
        # Plot the moment against time
        axs[0].plot(sol.t, moment, label='c = {}'.format(c))

    axs[0].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

TimeAndPositionOfModelOpening_plot_c(theta_i, theta_f, c_values_opening, k_start, k_finish)
TimeAndPositionOfModelClosing_plot_c(theta_f, theta_i, c_values_closing, k_start, k_finish)
