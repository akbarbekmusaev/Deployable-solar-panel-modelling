# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 15:17:35 2023

@author: hyper
"""

import numpy as np
from scipy.integrate import solve_ivp

def phase3(p, S2):
    #Step 1: Define initial conditions
    #S2y = [0.511468, 3.14159, -0.923012, 5.10871]
    #S2t = [1.59429]

    thet,phi,thetdot,phidot = S2.y
    t2 = S2.t
    
    g       = p[0]
    m_B     = p[1]
    M_CW    = p[2]
    M_P     = p[3]
    L_B     = p[4]
    H       = p[5]
    L_S     = p[6]
    L_BC    = p[7]
    L_BP    = p[4]-p[7]

    z2 = np.array([-L_BP * np.sin(thet[-1]) - L_S * np.sin(thet[-1] + np.pi - phi[-1]),
                   H + L_BP * np.cos(thet[-1]) + L_S * np.cos(thet[-1] + np.pi - phi[-1])])
    v2 = np.array([-L_BP * np.cos(thet[-1]) * thetdot[-1] - (L_S * np.cos(thet[-1] + np.pi - phi[-1]) * (thetdot[-1] - phidot[-1])),
                   -L_BP * np.sin(thet[-1]) * thetdot[-1] - (L_S * np.sin(thet[-1] + np.pi - phi[-1]) * (thetdot[-1] - phidot[-1]))])

    #testing
    #z2 = [-6.11822, 16.9003]
    #v2 = [36.7908, 20.6502]

    # Step 2: Define drag coefficient, density of air, and density of water
    Cd = 0.5  # drag coefficient
    rho_air = 1.225  # density of air
    rho_wat = 1000 # density of water

    # Step 3: Define drag force parameter DP, area of a sphere, radius

    r_p = ((3*M_P)/(4*rho_wat*np.pi))**(1/3)

    A_p = (np.pi)*(r_p)**2

    DP = 0.5 * Cd * rho_air * A_p  # Drag force parameter

    # Step 4: Define functions for drag forces in x- and y-directions
    def drag_force_x(vx, vy):
        # Define drag force in x-direction
        return DP * vx * (vx**2 + vy**2)**0.5

    def drag_force_y(vy, vx):
        # Define drag force in y-direction
        return DP * vy * (vx**2 + vy**2)**0.5

    # Step 5: Define the first-order equation of motion function f3
    def f3(t, z):
        x, y, vx, vy = z
        
        # Define drag forces in x- and y-directions
        F_drag_x = drag_force_x(vx, vy)
        F_drag_y = drag_force_y(vy, vx)
        
        # Define equations of motion
        dxdt = vx
        dydt = vy
        dvxdt = -M_P/(F_drag_x)
        dvydt = (-F_drag_y - M_P*g) / M_P # gravitational force
        
        return dxdt, dydt, dvxdt, dvydt

    # Step 6: Define the ODE solver options 'rtol', 'atol', and 'events'
    rtol = 1e-6  # Relative tolerance

    # Define the event function to stop integration when y-coordinate becomes zero
    def event_function(t, z):
        return z[1]

    event_function.terminal = True  # Stop integration when event occurs

    # Step 7: Set the initial conditions x, y, vx, vy
    initial_conditions = np.concatenate([z2, v2])
    # Append initial velocities to z2

    # Step 8: Run solve_ivp to find time (t3) and states (z3)
    S3 = solve_ivp(
        f3, [t2[-1],10000], initial_conditions,
        method='RK45', events = event_function, max_step = 0.1, rtol=rtol
    )

    # Extract results
    return S3

