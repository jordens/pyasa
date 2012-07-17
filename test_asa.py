# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Adaptive Simulated Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@gmail.com>

import numpy as np

import asa

s, t, c = .2, .05, .15
d = np.array([1., 1000., 10., 100.])

def cost(x):
    #raise ValueError("rrr")
    k = np.rint(x/s)
    r = np.fabs(k*s-x)
    p = np.sign(k)
    q = np.where(r<t, c*(p*p*k*s-p*t)**2, x**2)
    return (d*q).sum()

x0 = np.array([999., -1007, 1001, -903])
xmax = 1e4*np.ones_like(x0)

print asa.asa(cost, x0, -xmax, xmax, full_output=True)

# simple leak check
#while True:
#    try:
#        asa.asa(cost, x0, -xmax, xmax, full_output=True)
#    except Exception:
#        pass


d = np.arange(5.)
def cost(x):
    return -(np.sin(1./(x-d))*np.exp(-(x-d)**2)).sum()
x0 = -d
xmax = 1e1*np.ones_like(x0)

print asa.asa(cost, x0, -xmax, xmax, full_output=True)


