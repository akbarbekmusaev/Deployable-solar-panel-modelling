from Trebuchet_Functions import M1fun
from Trebuchet_Functions import R1fun, dv_Sfun
import numpy as np


def phase1(p):
    #input parameters
    g       = p[0]
    m_B     = p[1]
    M_CW    = p[2]
    M_P     = p[3]
    L_B     = p[4]
    H       = p[5]
    L_S     = p[6]
    L_BC    = p[7]
    L_BP    = p[8]


    #Setting up the initial thets
    thet_init = np.pi-np.arccos(H/L_BP)
    dthet_init = 0

    #Introduce solve_ivp
    from scipy.integrate import solve_ivp

    #the following function calculates ddthet from the inputs of thet and dthet
    def calculate_ddthet(thet, dthet):
        M = M1fun(thet, p)
        R = R1fun(thet, dthet, p)
        return -R/M

    #the following function is the function used by solve_ivp
    def next_step(t, thets):
        thet, dthet = thets
        ddthet = calculate_ddthet(thet, dthet)
        return (dthet, ddthet)

    #the following function is used with the event function to find out when Tv-mg=0
    def get_zero_Tvmg(t, thets):
       thet, dthet = thets
       ddthet = calculate_ddthet(thet, dthet)
       acceleration=dv_Sfun(thet, dthet, ddthet, p)
       Th=-M_CW * acceleration
       phi = np.arcsin((H-L_BP * np.sin(thet -np.pi/2))/L_S)+ thet - np.pi / 2
       Tv = Th * np.tan(np.pi/2 + phi - thets[0])
       Tvmg = Tv-M_CW * g
       return Tvmg
    get_zero_Tvmg.terminal = True

    #define thets
    thets = (thet_init, dthet_init) 

    #use solve_ivp to solve the ODE
    solution = solve_ivp(next_step, (0, 5), thets ,rtol=1e-6, events=get_zero_Tvmg)
    
    return solution
