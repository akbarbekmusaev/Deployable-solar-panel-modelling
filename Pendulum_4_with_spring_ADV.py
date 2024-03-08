# -*- coding: utf-8 -*-

# Setting initial angle to -pi/2 leads to error
# How many simulations are required to complete 20s?

import numpy as np
from scipy.integrate import solve_ivp
from Pendulum_Animations import animate_sprung_pendulum
import matplotlib.pyplot as plt

m = 1
g = 9.81
L = 0.75
k = 50
c = 1
Ls = 0.5
Ld = 0.6

def f(t,z):
    dz = [z[1],
          -g/L*np.sin(z[0])]
    return dz

def f_C(t,z):
    dz = [z[1],
          -g/L*np.sin(z[0]) - Ls**2*k/(m*L**2)*z[0] - Ld**2*c/(m*L**2)*z[1]]
    return dz

def my_event(t,z):
    return z[0]
my_event.terminal = True

# Solver params
T = 20
rtol = 1e-6
IC = [-np.pi/8, 0]

if IC[0] > 0:
    my_event.direction = -1
    sol = solve_ivp(f, (0,T), IC, rtol=rtol, events=my_event)
else:
    my_event.direction = +1
    sol = solve_ivp(f_C, (0,T), IC, rtol=rtol, events=my_event)

t_list = [sol.t]
y_list = [sol.y]

IC = sol.y[:,-1]
tspan = (t_list[-1][-1],T)

while t_list[-1][-1] < T:
    my_event.direction *= -1
    
    if my_event.direction == -1:
        sol = solve_ivp(f, tspan, IC, rtol=rtol, events=my_event)
    else:
        sol = solve_ivp(f_C, tspan, IC, rtol=rtol, events=my_event)
    
    t_list.append(sol.t)
    y_list.append(sol.y)
    
    IC = sol.y[:,-1]
    tspan = (t_list[-1][-1],T)

# # ANIMATION
# t = np.hstack(t_list)
# z = np.hstack(y_list)
# thet = z[0,:]

# ani = animate_sprung_pendulum(t,thet)

# PLOTTING
fig, ax1 = plt.subplots()

for i in range(len(t_list)):
    if abs(max(y_list[i][0])) > abs(min(y_list[i][0])):
        ax1.plot(t_list[i],y_list[i][0],'b')
    else:
        ax1.plot(t_list[i],y_list[i][0],'r')

ax1.set_xlabel(r'$t\;[s]$')
ax1.set_ylabel(r'$\theta\;[rad]$')

print(len(t_list))


