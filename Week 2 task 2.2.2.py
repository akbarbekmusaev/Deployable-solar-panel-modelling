import numpy as np
import matplotlib.pyplot as plt

#paramteres 
m = 600 #kg
P = 500 #N
c = 1 #Ns2m-2

#euler starting point
h = 0.1 #stepsize
v = 5 #ms-1
T = 15 #s
t = 0 #starting time in s

def p(t):
    p = (P/2)*(1+np.cos(np.pi*t))
    return p

def f(t, v): 
    f = (p(t) - c*v**2)/m
    return f

all_iters = []
all_iters.append(v)
data = np.zeros((2, int(T/h+1)))
data[0, :] = np.linspace(0, T, 151)

for i in range(0, int(T/h)):
    all_iters.append(all_iters[i] + h*f(t, all_iters[i]))
    t = t+h
    print(all_iters[i])


data[1, :] = all_iters
plt.plot(data[0, :], data[1, :])
maxvel = np.max(data, axis = 1)[1]

indices = np.where(data == maxvel)

print("max velocity = ")
print(maxvel)

print("max velocity time = ")
print(data[0, indices[1]])

