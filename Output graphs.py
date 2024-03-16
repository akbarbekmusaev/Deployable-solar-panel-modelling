from Initial_model import openingmodel
import matplotlib.pyplot as plt
import numpy as np

# Constants
theta_finishing = np.pi / 2
theta_initial = np.pi / 18
speed_initial = 0
T_stall = 0.43
omega_max = 2750
gear_ratios = [435,440,445,455]  # List of gear ratios

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
        sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratio)
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

# Main code
sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gear_ratios[0])  # Assuming gear_ratio[0] is used for simulation
TimeAndPositionOfModel_plot(sol, theta_finishing, gear_ratios)
total_powers = create_motor_plot(sol, gear_ratios)
plot_total_power(gear_ratios, total_powers)

