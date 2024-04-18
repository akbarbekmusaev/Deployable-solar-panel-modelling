from Initial_model import openingmodel
from Initial_model import closingmodel
from Initial_model import centre_of_mass
import matplotlib.pyplot as plt
import numpy as np

# Constants
theta_finishing = 100 * (np.pi / 180)
theta_initial = 5 * (np.pi / 180)
speed_initial = 0
T_stall = 0.69
omega_max = 3700
gear_ratios = [370]  # List of gear ratios
c = 0 # Damping coefficient
g = 9.81 # Acceleration due to gravity
M_total = 39.248 # Total mass of the pend]
efficiency = 0.73 # Efficiency of the gearbox



# Define function to plot time and position of model
def TimeAndPositionOfModel_plot(sol, theta_finishing, gear_ratios):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].set_title('Angular velocity of the model over Time')

    for gear_ratio in gear_ratios:
        sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio, c, 0, 0, efficiency)
        axs[0].plot(sol.t, sol.y[0, :], label='Gear Ratio {}'.format(gear_ratio))
        axs[1].plot(sol.t, sol.y[1, :], label='Gear Ratio {}'.format(gear_ratio))

    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].legend()
    axs[1].legend()
    plt.subplots_adjust(hspace=0.5)
    plt.show()

# Define function to plot motor values and return total power
def create_motor_plot(sol, gear_ratios):
    def torque_in(speed, gear_ratio):
        torque = -(T_stall / omega_max) * (speed * gear_ratio) + T_stall
        return torque
    
    total_powers = []
    for gear_ratio in gear_ratios:
        torque_values = [torque_in(speed, gear_ratio) for speed in sol.y[1, :]]
        motor_speed_values = [speed * gear_ratio for speed in sol.y[1, :]]
        power_values = [torque * (speed * (2 * np.pi / 60)) for torque, speed in zip(torque_values, motor_speed_values)]
        total_power = np.trapz(power_values, sol.t)
        print("Total power used during the operation (J) with gear ratio {}: {}".format(gear_ratio, total_power))
        total_powers.append(total_power)

    return total_powers

# Plot Total Power vs. Gear Ratio
def plot_total_power(gear_ratios, total_powers):
    plt.figure(figsize=(8, 6))
    plt.plot(gear_ratios, total_powers, marker='o')
    plt.xlabel('Gear Ratio')
    plt.ylabel('Total Power (J)')
    plt.title('Total Power used during the operation for different Gear Ratios')
    plt.grid(True)
    plt.show()
def holding_ratio(R_cm, theta):
    return (M_total * R_cm * np.cos(theta)*g)/T_stall
# Main code
sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratios[0], c, 0, 0, efficiency)  # Assuming gear_ratio[0] is used for simulation
TimeAndPositionOfModel_plot(sol, theta_finishing, gear_ratios)
# total_powers = create_motor_plot(sol, gear_ratios)
# plot_total_power(gear_ratios, total_powers)
# print("Centre of mass at initial angle: ", centre_of_mass(theta_initial))
# print("Holding ratio at initial angle: ", holding_ratio(centre_of_mass(theta_initial), theta_initial))

