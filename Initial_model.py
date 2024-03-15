import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum
from matplotlib.animation import PillowWriter
import webbrowser

def openingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio):
    # constants
    T = 10
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
        phi = theta + np.arcsin(L_base / L_shortbeam)
        return phi

    # Define function to calculate centre of mass
    def centre_of_mass(theta):
        # Calculate CoM for each component H is the horizontal distance from the pivot to the CoM
        # and V is the vertical distance

        # CoM of solar panels
        R_solar = np.sqrt((0.5 * L_solarpanel) ** 2 + L_longbeam ** 2 - L_solarpanel * L_longbeam * np.cos(phi(theta)))
        H_solar = L_longbeam * np.cos(theta) - 0.5 * L_solarpanel * np.cos(phi(theta) - theta)
        V_solar = np.sqrt(R_solar ** 2 - H_solar ** 2)
        angle_solar = np.arctan(V_solar / H_solar) * 180 / np.pi

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
            (0.5 * L_longbeam) ** 2 + (0.25 * L_longbeam) ** 2 - 0.25 * L_longbeam ** 2 * np.cos(np.pi - phi(theta)))
        H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
        V_green = np.sqrt(R_green ** 2 - H_green ** 2)
        angle_green = np.arctan(V_green / H_green) * 180 / np.pi

        # CoM of blue beam
        R_blue = np.sqrt(L_shortbeam ** 2 + (0.5 * L_longbeam) ** 2 - L_shortbeam ** 2 * np.cos(np.pi - phi(theta)))
        H_blue = np.sqrt(L_shortbeam ** 2 - L_base ** 2) + 0.5 * L_shortbeam * np.cos(theta)
        V_blue = np.sqrt(R_blue ** 2 - H_blue ** 2)
        angle_blue = np.arctan(V_blue / H_blue) * 180 / np.pi

        # H centre of mass
        H_cm = (
                           M_panels * H_solar + M_frame * H_frame + M_driving * H_driving + M_green * H_green + M_blue * H_blue) / M_total

        # V centre of mass
        V_cm = (
                           M_panels * V_solar + M_frame * V_frame + M_driving * V_driving + M_green * V_green + M_blue * V_blue) / M_total

        # R centre of mass
        R_cm = np.sqrt(H_cm ** 2 + V_cm ** 2)
        angle_cm = np.arctan(V_cm / H_cm) * 180 / np.pi

        return R_cm

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
    # Define function to calculate output torque from gearbox to mechanism
    def torque_out(speed):
        torque = torque_in(speed) * gearRatio
        return torque

    # Define function to calculate input torque from motor to gearbox
    def torque_in(speed):
        torque = -(T_stall / omega_max) * (speed * gearRatio) + T_stall
        return torque

    # Define function for differential equation
    def opening(t, z):
        difz = [z[1],
                (torque_out(z[1]) / (M_total * centre_of_mass(z[0]) ** 2)) - (g * np.cos(z[0])) / centre_of_mass(z[0])]
        return difz

    # Define event to stop the simulation when the pendulum reaches the finishing angle
    def my_eventstopopening(t, z):
        return z[0] - theta_finishing;

    my_eventstopopening.terminal = True
    my_eventstopopening.direction = 1

    # Solve the differential equation for opening
    sol_opening = solve_ivp(opening, (0, T), [theta_initial, speed_initial], rtol=0.000001, events=(my_eventstopopening))

    return sol_opening

def closingmodel(theta_initial, theta_finishing, speed_initial, T_stall, omega_max, gearRatio):
    # constants
    T = 10
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
        phi = theta + np.arcsin(L_base / L_shortbeam)
        return phi

    # Define function to calculate centre of mass
    def centre_of_mass(theta):
        # Calculate CoM for each component H is the horizontal distance from the pivot to the CoM
        # and V is the vertical distance

        # CoM of solar panels
        R_solar = np.sqrt((0.5 * L_solarpanel) ** 2 + L_longbeam ** 2 - L_solarpanel * L_longbeam * np.cos(phi(theta)))
        H_solar = L_longbeam * np.cos(theta) - 0.5 * L_solarpanel * np.cos(phi(theta) - theta)
        V_solar = np.sqrt(R_solar ** 2 - H_solar ** 2)
        angle_solar = np.arctan(V_solar / H_solar) * 180 / np.pi

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
            (0.5 * L_longbeam) ** 2 + (0.25 * L_longbeam) ** 2 - 0.25 * L_longbeam ** 2 * np.cos(np.pi - phi(theta)))
        H_green = 0.5 * L_longbeam * np.cos(theta) + 0.25 * L_longbeam * np.cos(phi(theta) - theta)
        V_green = np.sqrt(R_green ** 2 - H_green ** 2)
        angle_green = np.arctan(V_green / H_green) * 180 / np.pi

        # CoM of blue beam
        R_blue = np.sqrt(L_shortbeam ** 2 + (0.5 * L_longbeam) ** 2 - L_shortbeam ** 2 * np.cos(np.pi - phi(theta)))
        H_blue = np.sqrt(L_shortbeam ** 2 - L_base ** 2) + 0.5 * L_shortbeam * np.cos(theta)
        V_blue = np.sqrt(R_blue ** 2 - H_blue ** 2)
        angle_blue = np.arctan(V_blue / H_blue) * 180 / np.pi

        # H centre of mass
        H_cm = (
                           M_panels * H_solar + M_frame * H_frame + M_driving * H_driving + M_green * H_green + M_blue * H_blue) / M_total

        # V centre of mass
        V_cm = (
                           M_panels * V_solar + M_frame * V_frame + M_driving * V_driving + M_green * V_green + M_blue * V_blue) / M_total

        # R centre of mass
        R_cm = np.sqrt(H_cm ** 2 + V_cm ** 2)
        angle_cm = np.arctan(V_cm / H_cm) * 180 / np.pi

        return R_cm

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
                (-torque_out(z[1]) / (M_total * centre_of_mass(z[0]) ** 2)) - (g * np.cos(z[0])) / centre_of_mass(z[0])]
        return difz

    # Define event to stop the simulation when the pendulum reaches the finishing angle
    def my_eventstopclosing(t, z):
        return z[0] - theta_finishing;

    my_eventstopclosing.terminal = True
    my_eventstopclosing.direction = -1

    # Solve the differential equation for opening
    sol_closing = solve_ivp(closing, (0, T), [theta_initial, speed_initial], rtol=0.000001,events=(my_eventstopclosing))

    return sol_closing