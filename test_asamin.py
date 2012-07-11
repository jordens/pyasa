import numpy as np

from asa import asa

x0 = np.array([1.,])
xmin = np.array([0.,])
xmax = np.array([2.,])

def cost(x, xmin, xmax):
    print x, xmin, xmax
    return x**2, 0, 0

print asa(cost, x0, xmin, xmax)
