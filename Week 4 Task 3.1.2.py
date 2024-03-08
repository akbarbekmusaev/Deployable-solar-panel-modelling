import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from Pendulum_Animations import animate_pendulum


#constants
g = 9.81
l = 0.75
x0 = np.pi
v0 = 0
T = 2
iter_num = 100

def f(t, z):
    difz = [z[1], 
            -(g*np.sin(z[0]))/l]
    return difz

a = np.linspace(0.46, 0.47, iter_num)
finalVel = []
finalDis = []
cor_A = []

for i  in range(0, iter_num):
    sol = solve_ivp(f, (0, T), [a[i]*x0 ,v0], rtol = 0.00001)
    x = sol.y[0, :]
    v = sol.y[1, :]
    finalVel.append(v[-1])
    finalDis.append(x[-1])
    if round(v[-1], 3) == 0:
        cor_A.append(a[i])

plt.plot(a, finalVel)



#plt.plot(sol.t, x)

#ani = animate_pendulum(sol.t, x)

