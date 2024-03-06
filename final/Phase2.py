from Trebuchet_Functions import trebuchet_demo
import numpy as np
from scipy.integrate import solve_ivp
from Trebuchet_Functions import M2fun
from Trebuchet_Functions import R2fun
from Trebuchet_Functions import phifun
from Trebuchet_Functions import dphifun


def phase2(p, S1):
    p = [9.81, 12, 4000, 100, 8, 6, 6, 1.5, 6.5]

    t = S1.t[-1] #time from phase 1 ending
    thet_init = S1.y[0, -1]
    dthet_init = S1.y[1, -1]
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
    return solution

