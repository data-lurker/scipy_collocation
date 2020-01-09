import math as ma
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pylab as pl
import matplotlib.gridspec as gridspec
from scipy.optimize import fsolve

# x_dot
def f1(x):
    v = x[2]
    u = x[3]
    return v * np.sin(u)

# y_dot
def f2(x):
    v = x[2]
    u = x[3]
    return v * np.cos(u)

# v_dot
def f3(x):
    go = 9.81
    u = x[3]
    return go * np.cos(u)


brach_eom_list = [f1, f2, f3]


def objective_brach(x):
    return x[-1]


def ramp_initializer(x0, xf, y0, yf, v0, m, plot_validation):
    hk = 0.00001 # considered sufficient resolution for initialization fo this problem
    t = 0.0
    x = x0
    y = y0
    v = v0
    u = np.arctan((xf - x0) / (yf - y0))

    print "slope: {} radians | {} degrees".format(u, np.degrees(u))

    s = [x, y, v, u]

    xr = []
    yr = []
    vr = []
    tr = []

    xr.append(x)
    yr.append(y)
    vr.append(v)
    tr.append(t)

    c = 0

    while s[0] < xf:
        xo = xr[-1]
        yo = yr[-1]
        vo = vr[-1]

        s = [xo, yo, vo, u]

        x1 = xo + hk * f1(s)
        y1 = yo + hk * f2(s)
        v1 = vo + hk * f3(s)

        xr.append(x1)
        yr.append(y1)
        vr.append(v1)

        t += hk
        tr.append(t)

        c += 1
        
        if c>1e6:
            print("count exceeded 1e6")
            break 

    ti = np.linspace(0.0, t, m)
    ui = np.ones([m]) * u

    spl_x = InterpolatedUnivariateSpline(tr, xr)
    spl_y = InterpolatedUnivariateSpline(tr, yr)
    spl_v = InterpolatedUnivariateSpline(tr, vr)

    xi = spl_x(ti)
    yi = spl_y(ti)
    vi = spl_v(ti)

    xi[0] = x0
    xi[-1] = xf
    yi[0] = y0
    yi[-1] = yf
    vi[0] = v0

    xo = np.concatenate((xi, yi, vi, ui, [t]), axis=0).tolist()

    if plot_validation:
        gs = gridspec.GridSpec(2, 2)

        pl.figure()
        ax = pl.subplot(gs[0, 0])

        pl.plot(tr, xr, color="darkviolet")
        for i in range(len(ti)):
            pl.plot(ti[i], xi[i], marker="o", color="darkturquoise")
        pl.title("x position vs time")

        ax = pl.subplot(gs[0, 1])
        pl.plot(tr, yr, color="darkviolet")
        for i in range(len(ti)):
            pl.plot(ti[i], yi[i], marker="o", color="darkturquoise") 
        pl.title("y position vs time")           

        ax = pl.subplot(gs[1, 0])
        pl.plot(tr, vr, color="darkviolet")
        for i in range(len(ti)):
            pl.plot(ti[i], vi[i], marker="o", color="darkturquoise") 
        pl.title("v position vs time")                   

        ax = pl.subplot(gs[1, 1])
        pl.plot(xr, yr, color="darkviolet")
        for i in range(len(ti)):
            pl.plot(xi[i], yi[i], marker="o", color="darkturquoise")
        pl.title("x position vs y position")        

        pl.show()

    return xo


def slide_down_ramp(x0, xf, y0, yf, v0, hk):
    # hk = 0.001
    t = 0.0
    x = x0
    y = y0
    v = v0
    u = np.arctan((xf - x0) / (yf - y0))

    print "slope: {} radians | {} degrees".format(u, np.degrees(u))

    s = [x, y, v, u]

    xr = []
    yr = []
    vr = []

    xr.append(x)
    yr.append(y)
    vr.append(v)

    c = 0

    while s[0] < xf:
        xo = xr[-1]
        yo = yr[-1]
        vo = vr[-1]

        s = [xo, yo, vo, u]

        x1 = xo + hk * f1(s)
        y1 = yo + hk * f2(s)
        v1 = vo + hk * f3(s)

        xr.append(x1)
        yr.append(y1)
        vr.append(v1)

        t += hk
        c += 1
        
        if c>1e5:
            print("count exceeded 1e5")
            break 

    return (xr, yr, vr, t)


def cycloid_x(r, theta):
    return r * (theta - np.sin(theta))


def cycloid_y(r, theta):
    return r * (1.0 - np.cos(theta))


def cycloid_fxn_agr(p, xf, yf):
    r = p[0]
    theta = p[1]
    # xf = arg[0]
    # yf = arg[1]
    return [xf - cycloid_x(r, theta), yf - cycloid_y(r, theta)]


def cycloid_solver(xf, yf):
    p0 = [1.0, 2.0]
    r, theta = fsolve(func=cycloid_fxn_agr, x0=p0, args=(xf, yf))
    return (r, theta)


def cycloid_cordinates_xy(r, theta, n):
    xr = []
    yr = []

    for angle in np.linspace(0, theta, n):
        xr.append(cycloid_x(r, angle))
        yr.append(cycloid_y(r, angle))
    return (xr, yr)
    
    
def cycloid_time(xf, yf):
    r = cycloid_solver(xf, yf)[0]
    theta = cycloid_solver(xf, yf)[1]
    return theta * (r / 9.81) ** .5

