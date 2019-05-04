import numpy as np


def f1(x):
    v = x[4]
    u = x[6]
    return v * np.sin(u)


def f2(x):
    v = x[4]
    u = x[6]
    return v * np.cos(u)


def f3(x):
    go = 9.81
    u = x[6]
    return go * np.cos(u)


brach_eom_list = [f1, f2, f3]


def objective_brach(x):
    return x[-1]
