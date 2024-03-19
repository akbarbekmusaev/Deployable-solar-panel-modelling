from Initial_model import openingmodel
from Initial_model import closingmodel
import matplotlib.pyplot as plt
import numpy as np

# Constants
theta_finishing = np.pi / 18
theta_initial = np.pi / 2
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gear_ratio = 270
c = list(range(0, 100, 5))  # Damping coefficient


# Define function to plot time and position of model
def TimeAndPositionOfModel_plot(theta_finishing, c_values):
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

# Main code
# sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max,
#                    gear_ratio, c)
TimeAndPositionOfModel_plot(theta_finishing, c)
