# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 10:14:40 2023

@author: qq22099
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

g = 9.81 #acceleration due to gravity
Roh_A = 1.225 #density of air
Roh_W = 1000 #density of water
C_d = 0.5 #drag coefficient

p = [t2,z2]

#define the drag force parameter

#define functions for the drag forces in the x- and y- directions

#define the first order equation of motion, f3, as a funciton with inputs t and z using the drag force functions

# define the ODE solver options 'rol' and 'events'

#set the initial conditons x, y, xdot, ydot

#run solveivp to find time and states (t3,z3)

#no idea
def dd(x):
    return (-dr(x,y))/(M_p)
    
def dd(y):
    return -((M_p*g)-(dr(x,y)))/(M_p)
    
def dr(v):
    return D_P*v**2

def v(x,y):
    return (d(x)**2+d(y)**2)**0.5

def D_P(A_P):
    return 0.5*C_d*Roh_A*A_P
    
def A_P(r_P):
    return (np.pi)*(r_P)**2

def M_p(r_P):
    return Roh_W*(4/3)*(np.pi)*(r_P**3)
    