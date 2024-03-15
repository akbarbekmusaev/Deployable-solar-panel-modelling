from Initial_model import openingmodel
from Initial_model import closingmodel
import matplotlib.pyplot as plt
import numpy as np

#constants
theta_finishing = np.pi / 2
theta_initial = np.pi / 18
speed_initial = 0
T_stall = 0.19
omega_max = 5000
gearRatio = 1000

def TimeAndPositionOfModel_plot(sol, theta_finishing):
    fig, axs = plt.subplots(2, figsize=(10, 10))
    axs[0].plot(sol.t, sol.y[0, :], label='Angular displacement')
    axs[0].axhline(y=theta_finishing, color='r', linestyle='--', label='Finishing Theta')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Angular displacement (rad)')
    axs[0].legend()
    axs[0].set_title('Angular displacement of the model over Time')
    axs[1].plot(sol.t, sol.y[1, :], label='Angular velocity (rad/s)')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Angular velocity')
    axs[1].legend()
    axs[1].set_title('Angular velocity of the model over Time')
    plt.subplots_adjust(hspace=0.5)
    plt.show()

# Define function to plot motor values
def create_motor_plot(sol):
    def torque_in(speed):
        torque = -(T_stall / omega_max) * (speed * gearRatio) + T_stall
        return torque
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
    plt.subplots_adjust(hspace=1)
    plt.show()

sol = openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio)
create_motor_plot(sol)
TimeAndPositionOfModel_plot(sol, theta_finishing)
