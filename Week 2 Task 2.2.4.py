import numpy as np
from My_Euler_Solver import solve_Euler
import matplotlib.pyplot as plt

#import matplotlib.pyplot as plt
m = 600 #kg
P = 500 #N
c = 1 #Ns2m-2

def p(t, v):
    if v < 10:
        p = P
    else:
        p = (P/2)*(1+np.cos(np.pi*t))
    return p

def f(t, v): 
    f = (p(t, v) - c*v**2)/m
    return f

data = solve_Euler(f, 15, 5, 0.0001)
fig, ax1 = plt.subplots()
ax1.plot(data[0, :], data[1, :])
plt.xlim([6.5, 7.5])
plt.ylim([9.95, 10.05])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Velocity (m/s)')
ax1.set_title('Veloc vs time')
