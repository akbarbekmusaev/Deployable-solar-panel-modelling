import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

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

def spring_rate (theta, k_start, k_finish):
    # Define function to calculate the spring rate
    k = k_start + (k_finish - k_start) * (theta / np.pi/2)
    return k

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
    R_driving = 0.25 * L_longbeam
    H_driving = 0.25 * L_longbeam * np.cos(theta)
    V_driving = np.sqrt(R_driving ** 2 - H_driving ** 2)
    angle_driving = np.arctan(V_driving / H_driving) * 180 / np.pi

    # CoM of green beam
    R_green = np.sqrt(
        (0.5 * L_longbeam) ** 2 + (0.25 * L_longbeam) ** 2 - 0.25 * (L_longbeam ** 2) * np.cos(np.pi - phi(theta)))
    H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
    V_green = np.sqrt(R_green ** 2 - H_green ** 2)
    angle_green = np.arctan(V_green / H_green) * 180 / np.pi

    # Blue beam
    H_blue = L_shortbeam + 0.5 * L_shortbeam * np.cos(theta)
    V_blue = 0.5 * L_shortbeam * np.sin(theta) - L_base
    R_blue = np.sqrt(H_blue ** 2 + V_blue ** 2)
    angle_blue = np.arctan(V_blue / H_blue) * 180 / np.pi

    # H centre of mass
    H_cm = (M_solarpanel * (H_bottompanel + H_middlepanel + H_toppanel) + M_frame * H_frame + M_driving * H_driving + M_green * H_green + M_blue * H_blue) / M_total

    # V centre of mass
    V_cm = (M_solarpanel * (V_bottompanel + V_middlepanel + V_toppanel) + M_frame * V_frame + M_driving * V_driving + M_green * V_green + M_blue * V_blue) / M_total

    # R centre of mass
    R_cm = np.sqrt(H_cm ** 2 + V_cm ** 2)
    angle_cm = np.arctan(V_cm / H_cm) * 180 / np.pi

    return R_cm

def torque_out(speed, T_stall, omega_max, gearRatio):
    torque = torque_in(speed, T_stall, omega_max, gearRatio) * gearRatio
    return torque

def torque_in(speed, T_stall, omega_max, gearRatio):
    if speed < 0:
        speed = 0
    torque = -(T_stall / omega_max) * (speed * gearRatio) + T_stall
    return torque

def openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio, c, k_start, k_finish):
    # Define function for differential equation
    def opening(t, z):
        difz = [z[1],
                (torque_out(z[1], T_stall, omega_max, gearRatio) / (M_total * centre_of_mass(z[0]) ** 2)) - ((g * np.cos(z[0])) / centre_of_mass(z[0])) - c*z[1]/ (M_total * centre_of_mass(z[0]) ** 2) + spring_rate(z[0] - theta_initial, k_start, k_finish) / (M_total * centre_of_mass(z[0]) ** 2)]
        return difz

    # Define event to stop the simulation when the pendulum reaches the finishing angle
    def my_eventstopopening(t, z):
        return z[0] - theta_finishing;

    my_eventstopopening.terminal = True
    my_eventstopopening.direction = 1

    # Solve the differential equation for opening
    sol_opening = solve_ivp(opening, (0, T), [theta_initial, speed_initial], rtol=0.000001, events=my_eventstopopening)

    return sol_opening


def closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio, c, k_start, k_finish):
    # Define function for differential equation
    def closing(t, z):
        difz = [z[1],
                (-torque_out(z[1], T_stall, omega_max, gearRatio) / (M_total * centre_of_mass(z[0]) ** 2)) - (g * np.cos(z[0])) / centre_of_mass(z[0])  - c*z[1]/ (M_total * centre_of_mass(z[0]) ** 2) + spring_rate(theta_finishing - z[0], k_finish, k_start) / (M_total * centre_of_mass(z[0]) ** 2)]
        return difz

    # Define event to stop the simulation when the pendulum reaches the finishing angle
    def my_eventstopclosing(t, z):
        return z[0] - theta_finishing;

    my_eventstopclosing.terminal = True
    my_eventstopclosing.direction = -1

    # Solve the differential equation for opening
    sol_closing = solve_ivp(closing, (0, T), [theta_initial, speed_initial], rtol=1e-6, events=my_eventstopclosing)

    return sol_closing
