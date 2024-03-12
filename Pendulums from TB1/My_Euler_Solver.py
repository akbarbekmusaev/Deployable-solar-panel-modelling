import numpy as np

def solve_Euler(f, T, v0, h):

    t = 0 #starting time in s
    
    all_iters = []
    all_iters.append(v0)
    data = np.zeros((2, int(T/h+1)))
    data[0, :] = np.linspace(0, T, int(T/h+1))

    for i in range(0, int(T/h)):
        all_iters.append(all_iters[i] + h*f(t, all_iters[i]))
        t = t+h
        
    data[1, :] = all_iters
    return data