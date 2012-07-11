# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Adaptive Simulated Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@gmail.com>

import numpy as np

import asa

def asa_test_cost(x):
    s_i = .2
    t_i = .05
    c_r = .15
    q_n = 0.
    for i in range(x.shape[0]):
        d_i = [1., 1000., 10., 100][i]
        if x[i] > 0:
            k_i = int(x[i] / s_i + .5)
        elif x[i] < 0:
            k_i = int(x[i] / s_i - .5)
        else:
            k_i = 0
        if abs(k_i * s_i - x[i]) < t_i:
            if k_i < 0:
                z_i = k_i * s_i + t_i
            elif k_i > 0:
                z_i = k_i * s_i - t_i
            else:
                z_i = 0.
            q_n += c_r * d_i * z_i * z_i
        else:
            q_n += d_i * x[i] * x[i]
    return q_n

x0 = np.array([999., -1007, 1001, -903])
n = x0.shape[0]
xmax = 1e4*np.ones((4,))
xmin = -xmax

print asa.asa(asa_test_cost, x0, xmin, xmax, full_output=True)
