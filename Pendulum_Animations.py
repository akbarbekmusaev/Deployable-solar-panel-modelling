# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d

##########################################################################################

def animate_pendulum(t,thet):
    # Default FPS
    FPS = 30
    
    # Remove non-unique instances
    t,idx = np.unique(t,return_index=True)
    thet = thet[idx]
    
    # Relinearise time and displacement interpolate
    thetf = interp1d(t, thet, kind='cubic')
    tL = np.linspace(t[0],t[-1],round((t[-1]-t[0])*FPS+1))
    thetL = thetf(tL)
    
    L = 0.75
    x = L*np.sin(thetL)
    y = -L*np.cos(thetL)
    
    # Animation function
    def animate(i):
        plh_arm.set_data((0,x[i]), (0,y[i]))
        plh_bob.set_data(x[i], y[i])
    
    # Initialise plot
    fig, ax = plt.subplots()
    plh_arm, = ax.plot([], [], color='blue', linewidth=2.0)
    plh_bob, = ax.plot([], [], marker='.', color='blue', markersize=20)
    ax.plot(0,0,marker='.',color='red',markersize=8)
    
    ax.set_xlim([-1,+1])
    ax.set_ylim([-1,+1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    
    # Create animation
    int_time = 1000*(tL[-1]-tL[0])/(len(tL)-1)     # Interval time [ms]
    ani = FuncAnimation(fig=fig,
                        func=animate, 
                        frames=range(len(tL)), 
                        interval=int_time, 
                        repeat=False)
    return ani

###############################################################################

def animate_sprung_pendulum(t,thet):
    # Default FPS
    FPS = 30
    
    # Remove non-unique instances
    t,idx = np.unique(t,return_index=True)
    thet = thet[idx]
    
    # Relinearise time and displacement interpolate
    thetf = interp1d(t, thet, kind='cubic')
    tL = np.linspace(t[0],t[-1],round((t[-1]-t[0])*FPS+1))
    thetL = thetf(tL)
    
    L = 0.75
    x = L*np.sin(thetL)
    y = -L*np.cos(thetL)
    
    # Create masks for spring positions
    N = 10
    spr_thet0 = np.linspace(-np.pi,-7*np.pi/8,100)
    spr_thet1M = np.linspace(0,1,4*N+1)
    spr_thet1M = spr_thet1M[1::2]
    spr_thetNM = np.linspace(0,1,100)
    
    spr_r0 = 0.5*np.ones(spr_thet0.shape)
    spr_r1 = np.resize([0.55,0.45], len(spr_thet1M))
    spr_rN = 0.5*np.ones(spr_thetNM.shape)
    spr_r = np.concatenate((spr_r0,spr_r1,spr_rN,[0.45,0.55]))
    
    # Animation function
    def animate(i):
        plh_arm.set_data((0,x[i]), (0,y[i]))
        plh_bob.set_data(x[i], y[i])
        
        # Compute spring positions
        if thetL[i] < 0:
            spr_T0 = thetL[i]
        else:
            spr_T0 = 0
        
        spr_T1 = spr_T0-1*np.pi/8
        
        spr_thet1 = (-7*np.pi/8) + spr_thet1M*( spr_T1 - (-7*np.pi/8) )
        spr_thetN = (spr_T1) + spr_thetNM*( spr_T0 - spr_T1 )
        
        spr_thet = np.concatenate((spr_thet0,spr_thet1,spr_thetN,[spr_T0,spr_T0]))
        
        spr_x = +spr_r*np.sin(spr_thet)
        spr_y = -spr_r*np.cos(spr_thet)
        
        plh_spr.set_data(spr_x,spr_y)
    
    # Initialise plot
    fig, ax = plt.subplots()
    
    # Create pendulum features
    plh_arm, = ax.plot([], [], color='blue', linewidth=2.0)
    plh_bob, = ax.plot([], [], marker='.', color='blue', markersize=20)
    ax.plot(0,0,marker='.',color='red',markersize=8)
    
    # Create spring features
    plh_spr, = ax.plot([],[],color='blue')
    ax.plot([0,0],[0.45,0.55],'b',linewidth=4)
    
    ax.set_xlim([-1,+1])
    ax.set_ylim([-1,+1])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    
    # Create animation
    int_time = 1000*(tL[-1]-tL[0])/(len(tL)-1)     # Interval time [ms]
    ani = FuncAnimation(fig=fig,
                        func=animate, 
                        frames=range(len(tL)), 
                        interval=int_time, 
                        repeat=False)
    
    return ani