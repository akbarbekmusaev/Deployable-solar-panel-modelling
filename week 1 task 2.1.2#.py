import numpy as npfd
import time


h = [10**-1, 10**-2, 10**-3, 10**-4, 1.000001*10**-5] #stepsize
fig, ax1 = plt.subplots()

final_vel = []
log_h = []
exec_time = []

for i in range(0, len(h)):
    tS = time.perf_counter()
    #paramteres 
    m = 600 #kg
    P = 500 #N
    c = 1 #Ns2m-2

    #euler starting point
    v = 5 #ms-1
    T = 15 #s
    t = 0 #starting time in s
    
    all_iters = []
    data = np.zeros((2, int(T/h[i]+1)))
    data[0, :] = np.linspace(0, 15, int(T/h[i]+1))
    while t <= T:
        v = v + (h[i]/m)*((P/2)*(1+np.cos(np.pi*t))-c*v**2)
        all_iters.append(v)
        t = t + h[i]

    data[1, :] = all_iters
    
    final_vel.append(data[1, int(T/h[i])])
    log_h.append(np.log(h[i]))
    
    ax1.plot(data[0, :], data[1, :], label = h[i])
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Velocity (m/s)')
    ax1.set_title('Veloc vs time')
    ax1.legend()

    maxvel = np.max(data, axis = 1)[1]
    indices = np.where(data == maxvel)
    
    print("max velocity = " + str(maxvel))
    print("max velocity time = " + str(data[0, indices[1]]))
    tE = time.perf_counter()
    exec_time.append(tE-tS)


fig, ax2 = plt.subplots()
ax2.plot(log_h, final_vel)
ax2.set_title('fin vel vs log h')

fig, ax3 = plt.subplots()
ax3.plot(log_h, exec_time)
ax3.set_title('exec time vs log h')