# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 10:00:19 2023

@author: hyper
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Initial conditions
t2 = 0
z2 = [0, 2]

# Step 1: Define parameters
def define_parameters(p):
    return [t2, z2]

# Step 2: Define drag coefficient, density of air, and density of water
Cd = 0.5  # drag coefficient
rho_air = 1.225  # density of air

# Step 3: Define drag force parameter DP
DP = 0.5 * Cd * rho_air  # Drag force parameter

# Step 4: Define functions for drag forces in x- and y-directions
def drag_force_x(vx, vy):
    # Define drag force in x-direction
    return -DP * vx * (vx**2 + vy**2)**0.5

def drag_force_y(vy, vx):
    # Define drag force in y-direction
    return -DP * vy * (vx**2 + vy**2)**0.5

# Step 5: Define the first-order equation of motion function f3
def f3(t, z):
    x, y, vx, vy = z
    
    # Define drag forces in x- and y-directions
    F_drag_x = drag_force_x(vx, vy)
    F_drag_y = drag_force_y(vy, vx)
    
    # Define equations of motion
    dxdt = vx
    dydt = vy
    dvxdt = -F_drag_x
    dvydt = -F_drag_y - 9.81  # gravitational force
    
    return [dxdt, dydt, dvxdt, dvydt]

# Step 6: Define the ODE solver options 'rtol', 'atol', and 'events'
rtol = 1e-6  # Relative tolerance

ode_solver_options = {'rtol': rtol}

# Step 7: Set the initial conditions x, y, vx, vy
initial_conditions = z2 + [2, 1]  # Append initial velocities to z2

# Step 8: Run solve_ivp to find time (t3) and states (z3)
solution = solve_ivp(
    f3, [t2, np.inf], initial_conditions,
    method='RK45', dense_output=True, **ode_solver_options
)

# Extract results
t3 = solution.t
z3 = solution.y

# Function to update the plot in each animation frame
def update(frame):
    plt.clf()
    
    # Plot the trajectory up to the current frame
    plt.plot(z3[0, :frame], z3[1, :frame], 'b-', label='Projectile Motion')
    
    # Plot the current position of the projectile
    plt.plot(z3[0, frame], z3[1, frame], 'ro', label='Current Position')
    
    # Draw the floor (horizontal line at y=0)
    plt.axhline(y=0, color='k', linestyle='--', label='Ground')

    plt.title(f'Projectile Motion at t = {t3[frame]:.2f} seconds')
    plt.xlabel('Horizontal Distance (meters)')
    plt.ylabel('Vertical Distance (meters)')
    plt.legend()
    
    # Print velocities at the current frame
    print(f"At t = {t3[frame]:.2f} seconds: vx = {z3[2, frame]:.4f}, vy = {z3[3, frame]:.4f}")

# Create a figure and axis for the animation
fig, ax = plt.subplots()

# Set the axis limits based on the maximum distances in the x and y directions
ax.set_xlim(0, np.max(z3[0, :]))
ax.set_ylim(-10, np.max(z3[1, :]))  # Adjust the y-axis limits as needed

# Create the animation
animation = FuncAnimation(fig, update, frames=len(t3), interval=50, repeat=False)

# Show the animation
plt.show()