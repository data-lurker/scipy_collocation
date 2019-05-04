import math as ma
import numpy as np
from scipy.optimize import minimize
from problem_functions import brach_eom
from utilities import collocation_functions as cf

ns = 5
m = ns + 1
p = 6
q = 1
gamma = p + q

# initialize
t1 = 0.0
tf = .4
hk = (tf - t1) / ns
t = np.linspace(t1, tf, 100)

# Boundry Conditions / Constrains
# Initial
x0 = 0.0
xd0 = 0.0
y0 = 0.0
yd0 = 0.0
v0 = 0.0
vd0 = 0.0
# t0 = 0.0

# Final
xf = 1.0
yf = -1.0
tf = 3.0

x = np.linspace(x0, xf, m)
xd = np.linspace(xd0, 10, m)
y = np.linspace(y0, yf, m)
yd = np.linspace(yd0, 10, m)
v = np.linspace(v0, 10, m)
vd = np.linspace(yd0, 10, m)
u = np.ones([m]) * np.tan((yf - y0) / (xf - x0))
# t = np.linspace(t0, tf, m)

# np.random.uniform(-.05, .05, m)

xo = np.concatenate((x, xd, y, yd, v, vd, u, [tf]), axis=0)

print(xo)
print(len(xo))
