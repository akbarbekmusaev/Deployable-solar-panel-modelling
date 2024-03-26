import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import simpy as sp

# constants
T = 100
g = 9.81  # m/s^2

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
    return theta + np.arcsin(L_base / L_shortbeam)

def damperposition(theta):
    # End-point
    H = 0.5 * L_longbeam * np.cos(theta) + 0.5 * L_longbeam * np.cos(phi(theta) - theta)
    V = 0.5 * L_longbeam * np.sin(theta) - 0.5 * L_longbeam * np.sin(phi(theta) - theta) + L_base
    R = np.sqrt(H ** 2 + V ** 2)
    # Blue beam
    R_blue = np.sqrt(L_shortbeam ** 2 + (0.5 * L_shortbeam) ** 2 - L_shortbeam ** 2 * np.cos(np.pi - phi(theta)))
    H_blue = np.sqrt(L_shortbeam ** 2 - L_base ** 2) + 0.5 * L_shortbeam * np.cos(theta)
    V_blue = np.sqrt(R_blue ** 2 - H_blue ** 2)
    # Green beam
    R_green = np.sqrt(
        (0.5 * L_longbeam) ** 2 + (0.25 * L_longbeam) ** 2 - 0.25 * (L_longbeam ** 2) * np.cos(np.pi - phi(theta)))
    H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
    V_green = np.sqrt(R_green ** 2 - H_green ** 2)
    # Angle a
    H_a = H - H_blue
    V_a = V - V_blue
    R_a = np.sqrt(H_a ** 2 + V_a ** 2)
    angle_a = np.arccos(((H_a) * H + (V_a) * V) / (R_a * R))
    # angle b
    H_b = H - H_green
    V_b = V - V_green
    R_b = np.sqrt(H_b ** 2 + V_b ** 2)
    angle_b = np.arccos(((H_b) * H + (V_b) * V) / (R_b * R))
    return H, V, H_green, V_green, H_blue, V_blue

def dampertorque(c, theta, omega):
    omega = abs(omega)
    angle_a, angle_b = damperposition(theta)
    #Calculate speed and torque
    H_prime = -0.5 * omega * L_longbeam * np.sin(theta)
    V_prime = 0.5 * omega * L_longbeam * np.cos(theta)
    speed = np.sqrt(H_prime ** 2 + V_prime ** 2)
    F_magnitude = c * speed
    F_green = F_magnitude * (np.sin(angle_b)/np.sin(angle_a + angle_b))
    Torque = F_green * L_shortbeam * abs(np.sin(phi(theta)))

    return Torque

def springtorque(k, theta, theta_initial):
    #initial position of the spring
    H_initial = 0.5 * L_longbeam * np.cos(theta_initial) + 0.5 * L_longbeam * np.cos(
        phi(theta_initial) - theta_initial)
    V_initial = 0.5 * L_longbeam * np.sin(theta_initial) - 0.5 * L_longbeam * np.sin(
        phi(theta_initial) - theta_initial)
    R_initial = np.sqrt(H_initial ** 2 + V_initial ** 2)
    # End-point
    H = 0.5 * L_longbeam * np.cos(theta) + 0.5 * L_longbeam * np.cos(phi(theta) - theta)
    V = 0.5 * L_longbeam * np.sin(theta) - 0.5 * L_longbeam * np.sin(phi(theta) - theta)
    R = np.sqrt(H ** 2 + V ** 2)
    angle_a, angle_b = damperposition(theta)

    # Calculate displacement and torque
    displacement = abs(R_initial - R)
    F_magnitude = k * displacement

    F_green = F_magnitude * (np.sin(angle_b) / np.sin(angle_a + angle_b))
    Torque = F_green * L_shortbeam * np.cos(phi(theta) - np.pi / 2)

    return Torque

def springtorqueconstant(k, theta, theta_initial):
    angle_a, angle_b = damperposition(theta)

    # Calculate displacement and torque
    F_magnitude = k

    F_green = F_magnitude * (np.sin(angle_b) / np.sin(angle_a + angle_b))
    Torque = F_green * L_shortbeam * np.cos(abs(phi(theta)))

    return Torque


def centre_of_mass(theta):
    # Calculate CoM for each component H is the horizontal distance from the pivot to the CoM
    # and V is the vertical distance

    # CoM of solar panels
    H_toppanel =1.5 * L_longbeam * np.cos(theta) - 0.5 * L_longbeam * np.cos(phi(theta) - theta)-0.5*L_solarpanel*np.cos(phi(theta)-theta)
    V_toppanel = 1.5 * L_longbeam * np.sin(theta) + 0.5 * L_longbeam * np.sin(phi(theta) - theta)+0.5*L_solarpanel*np.sin(phi(theta)-theta)
    R_toppanel = np.sqrt(H_toppanel**2+V_toppanel**2)
    angle_toppanel = np.arctan(V_toppanel / H_toppanel) * 180 / np.pi

    H_middlepanel = L_longbeam * np.cos(theta) - 0.5 * L_solarpanel * np.cos(phi(theta) - theta)
    V_middlepanel = L_longbeam * np.sin(theta) + 0.5 * L_solarpanel * np.sin(phi(theta) - theta)
    R_middlepanel = np.sqrt(H_middlepanel**2+V_middlepanel**2)
    angle_middlepanel = np.arctan(V_middlepanel / H_middlepanel) * 180 / np.pi

    H_bottompanel = 0.5 * L_longbeam * np.cos(theta) + (0.5 * L_longbeam - 0.5*L_solarpanel)*np.cos(phi(theta)-theta)
    V_bottompanel = 0.5 * L_longbeam * np.sin(theta) - (0.5 * L_longbeam - 0.5*L_solarpanel)*np.sin(phi(theta)-theta)
    R_bottompanel = np.sqrt(H_bottompanel**2+V_bottompanel**2)
    angle_bottompanel = np.arctan(V_bottompanel / H_bottompanel) * 180 / np.pi

    # CoM of top frame
    R_frame = np.sqrt((0.5 * L_longbeam) ** 2 + L_longbeam ** 2 - L_longbeam ** 2 * np.cos(phi(theta)))
    H_frame = L_longbeam * np.cos(theta) - 0.5 * L_longbeam * np.cos(phi(theta) - theta)
    V_frame = np.sqrt(R_frame ** 2 - H_frame ** 2)
    angle_frame = np.arctan(V_frame / H_frame) * 180 / np.pi

    # CoM of driving beam
    R_driving = 0.25 * L_longbeam * np.sin(theta) + L_base
    H_driving = 0.25 * L_longbeam * np.cos(theta)
    V_driving = np.sqrt(R_driving ** 2 - H_driving ** 2)
    angle_driving = np.arctan(V_driving / H_driving) * 180 / np.pi

    # CoM of green beam
    R_green = np.sqrt(
        (0.5 * L_longbeam) ** 2 + (0.25 * L_longbeam) ** 2 - 0.25 * (L_longbeam ** 2) * np.cos(np.pi - phi(theta)))
    H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
    V_green = np.sqrt(R_green ** 2 - H_green ** 2)
    angle_green = np.arctan(V_green / H_green) * 180 / np.pi

    # CoM of blue beam
    H_blue = np.sqrt(L_shortbeam ** 2 - L_base ** 2) + L_shortbeam * np.cos(theta)
    V_blue = 0.5 * np.sin(phi(theta)) * L_shortbeam
    R_blue = np.sqrt(H_blue ** 2 + V_blue ** 2)
    angle_blue = np.arctan(V_blue / H_blue) * 180 / np.pi

    # H centre of mass
    H_cm = (M_solarpanel * (H_toppanel + H_middlepanel + H_bottompanel) + M_frame * H_frame + M_driving * H_driving + M_green * H_green + M_blue * H_blue) / M_total

    # V centre of mass
    V_cm = (M_solarpanel * (V_toppanel + V_middlepanel + V_bottompanel) + M_frame * V_frame + M_driving * V_driving + M_green * V_green + M_blue * V_blue) / M_total

    # R centre of mass
    R_cm = np.sqrt(H_cm ** 2 + V_cm ** 2)
    angle_cm = np.arctan(V_cm / H_cm) * 180 / np.pi

    return H_blue, V_blue, H_green, V_green, H_driving, V_driving, H_frame, V_frame, H_toppanel, V_toppanel, H_middlepanel, V_middlepanel, H_bottompanel, V_bottompanel


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


def openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio, c, k):

    # Define function to calculate output torque from gearbox to mechanism
    def torque_out(speed):
        torque = torque_in(speed) * gearRatio
        return torque

    # Define function to calculate input torque from motor to gearbox
    def torque_in(speed):
        if speed < 0:
            speed = 0
        torque = -(T_stall / omega_max) * (speed * gearRatio) + T_stall
        return torque

    # Define function for differential equation
    def opening(t, z):
        difz = [z[1],
                (torque_out(z[1]) / (M_total * centre_of_mass(z[0]) ** 2)) - ((g * np.cos(z[0])) / centre_of_mass(z[0])) - (dampertorque(c, z[0], z[1])/(M_total*centre_of_mass(z[0]) ** 2))]
        return difz

    # Define event to stop the simulation when the pendulum reaches the finishing angle
    def my_eventstopopening(t, z):
        return z[0] - theta_finishing;

    my_eventstopopening.terminal = True
    my_eventstopopening.direction = 1

    # Solve the differential equation for opening
    sol_opening = solve_ivp(opening, (0, T), [theta_initial, speed_initial], rtol=0.000001, events=my_eventstopopening)

    return sol_opening


def closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio, c, k):
    # Define function to calculate output torque from gearbox to mechanism
    def torque_out(speed):
        torque = torque_in(speed) * gearRatio
        return torque

    # Define function to calculate input torque from motor to gearbox
    def torque_in(speed):
        torque = -(T_stall / omega_max) * (speed * gearRatio) + T_stall
        return torque

    # Define function for differential equation
    def closing(t, z):
        difz = [z[1],
                (-torque_out(z[1]) / (M_total * centre_of_mass(z[0]) ** 2)) - (g * np.cos(z[0])) / centre_of_mass(z[0]) + (dampertorque(c, z[0], z[1]) / (M_total*centre_of_mass(z[0]) ** 2)) + (springtorqueconstant(k, z[0], theta_finishing) / (M_total*centre_of_mass(z[0]) ** 2))]
        return difz

    # Define event to stop the simulation when the pendulum reaches the finishing angle
    def my_eventstopclosing(t, z):
        return z[0] - theta_finishing;

    my_eventstopclosing.terminal = True
    my_eventstopclosing.direction = -1

    # Solve the differential equation for opening
    sol_closing = solve_ivp(closing, (0, T), [theta_initial, speed_initial], rtol=1e-6, events=my_eventstopclosing)

    return sol_closing

# # Convert theta range from degrees to radians
# theta_values = np.radians(np.arange(15, 100))
#
# # Initialize lists to store H and V values
# H_values = [[] for _ in range(7)]
# V_values = [[] for _ in range(7)]
#
# # Calculate H and V values for each theta
# for theta in theta_values:
#     H_V_values = centre_of_mass(theta)
#     for i in range(7):
#         H_values[i].append(H_V_values[i*2])
#         V_values[i].append(H_V_values[i*2+1])
#
# # Plot H and V values
# plt.figure(figsize=(10, 6))
# labels = ['Blue', 'Green', 'Driving', 'Frame', 'Top Panel', 'Middle Panel', 'Bottom Panel']
# for i in range(7):
#     plt.plot(H_values[i], V_values[i], label=labels[i])
# plt.xlabel('H')
# plt.ylabel('V')
# plt.legend()
# plt.title('H vs V for different components and Centre of Mass')
# plt.show()

# Convert theta range from degrees to radians
theta_values = np.radians(np.arange(15, 101))

# Initialize lists to store H and V values
H_regular, V_regular, H_blue, V_blue, H_green, V_green = [], [], [], [], [], []

# Calculate H and V values for each theta
for theta in theta_values:
    H_reg, V_reg, H_gre, V_gre, H_blu, V_blu = damperposition(theta)
    H_regular.append(H_reg)
    V_regular.append(V_reg)
    H_blue.append(H_blu)
    V_blue.append(V_blu)
    H_green.append(H_gre)
    V_green.append(V_gre)

# Plot H and V values
plt.figure(figsize=(10, 6))
plt.plot(H_regular, V_regular, label='Point')
plt.plot(H_blue, V_blue, label='Blue')
plt.plot(H_green, V_green, label='Green')
plt.xlabel('H')
plt.ylabel('V')
plt.legend()
plt.title('H vs V for Regular, Blue, and Green')
plt.show()

