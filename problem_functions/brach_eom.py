import numpy as np


def f1(x):
    v = x[2]
    u = x[3]
    return v * np.sin(u)


def f2(x):
    v = x[2]
    u = x[3]
    return v * np.cos(u)


def f3(x):
    go = 9.81
    u = x[3]
    return go * np.cos(u)


brach_eom_list = [f1, f2, f3]


def objective_brach(x):
    return x[-1]


def initializer(x0, xf, y0, yf, v0, ns):
    hk = 0.001
    t = 0.0
    x = x0
    y = y0
    v = v0
    u = np.tan((yf - y0) / (xf - x0))

    s = [x, y, v, u]

    xr = []
    yr = []
    vr = []

    while s[0] < xf:
        s[0] += s[0] + hk * f1(s)
        s[1] += s[1] + hk * f2(s)
        s[2] += s[2] + hk * f3(s)
        xr.append(s[0])
        yr.append(s[1])
        vr.append(s[2])
        t += hk

    print("iteration compleate. tf={}".format(t))

    


    
