from Trebuchet_Functions import trebuchet_demo
import numpy as np
from scipy.integrate import solve_ivp
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


S1, S2, S3 = trebuchet_demo(H, L_S, L_BC)

t = 0.72517467 #time from phase 1 ending
thet_init = 2.20412634 
dthet_init = -1.73017
dphi_init = dphifun(thet_init, dthet_init, p)
phi_init = phifun(thet_init, p)

def f2(t,z):
    thet,phi,dthet,dphi=z
    M2 = M2fun(phi, p)
    R2 = R2fun(thet, phi, dthet, dphi, p)
    (ddthet,ddphi)=-((np.linalg.inv(M2)) @ (R2))  
    return (dthet,dphi,ddthet, ddphi)

def get_zero_phipi(t,z):
    phi_pi=z[1]-np.pi
    return phi_pi
get_zero_phipi.terminal = True

z=(thet_init,phi_init,dthet_init,dphi_init)

solution = solve_ivp(f2, (t, 10), z ,rtol=1e-6,events= get_zero_phipi)