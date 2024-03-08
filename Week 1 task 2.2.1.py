import numpy as np
import matplotlib.pyplot as plt


def y(a, b, c, x):
    return a*x**2 + b*x + c;

x = np.array(range(11))
a = 2
b = -20
c = 4
plt.plot(x, y(a, b, c, x))
plt.plot(x, y(1, -10, 4, x))