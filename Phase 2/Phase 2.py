from Trebuchet_Functions import trebuchet_demo
from Trebuchet_Functions import animate_trebuchet

p = [9.81, 12, 4000, 100, 8, 6, 6, 1.5, 6.5]

g       = p[0]
m_B     = p[1]
M_CW    = p[2]
M_P     = p[3]
L_B     = p[4]
H       = p[5]
L_S     = p[6]
L_BC    = p[7]
L_BP    = p[8]


S1, S2, S3 = trebuchet_demo(H, L_S, L_BC)


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Trebuchet_Functions import M2fun
from Trebuchet_Functions import R2fun
from Trebuchet_Functions import phifun
from Trebuchet_Functions import dphifun

p = [9.81, 12, 4000, 100, 8, 6, 6, 1.5, 6.5]

g       = p[0]
m_B     = p[1]
M_CW    = p[2]
M_P     = p[3]
L_B     = p[4]
H       = p[5]
L_S     = p[6]
L_BC    = p[7]
L_BP    = p[8]

t =  S1.t
thet =S1.y[0, 7]
dthet = S1.y[1,7]

rtol = 1e-6

def M2(z):
    M2[0] = M2fun(phi, p)
    return M2[0]

def R2(z):
    R2[0] = R2fun(thet, phi, dthet, dphi, p)
    return R2[0]

def phi(z):
    phi[0] = phifun(thet, p)
    return phi[0]

def dphi(z):
    dphi[0] = dphifun(thet, dthet, p)
    return dphi[0]
def f2(t, z):
    thet, phi, dthet, dphi = z
    M = M2(z)
    R = R2(z)
    dthet = -R/M
    return [thet, phi, dthet, dphi]


def P2event(t, z):
    phi = z
    return phi - np.pi
P2event.terminal = True
P2event.direction = 1

T = 100
z0 = [thet, dthet, phi, dphi]

sol2 = solve_ivp(f2, (0,T), z0, rtol=rtol, events=P2event)
t2 = sol2.t
z2 = sol2.y

    
