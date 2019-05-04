import math as ma
import numpy as np
from scipy.optimize import minimize
from problem_functions import brach_eom as beom
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

# Boundry Conditions / Constraints
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

tf0 = 3.0

x = np.linspace(x0, xf, m)
xd = np.linspace(xd0, 10, m)
y = np.linspace(y0, yf, m)
yd = np.linspace(yd0, 10, m)
v = np.linspace(v0, -10, m)
vd = np.linspace(yd0, 10, m)
u = np.ones([m]) * np.tan((yf - y0) / (xf - x0))
# np.random.uniform(-.05, .05, m)

xo = np.concatenate((x, xd, y, yd, v, vd, u, [tf0]), axis=0).tolist()

print(xo)
print(len(xo))

collocation_cons = cf.constraint_list_gen(ns, p, q, cf.dfct_eval, beom.brach_eom_list)

print(len(collocation_cons))
for i in range(len(collocation_cons)):
    print(collocation_cons[i])


# constraints
def c_x0(x):
    return x0 - x[0 * m]


def c_xd0(x):
    return xd0 - x[1 * m]


def c_y0(x):
    return y0 - x[2 * m]


def c_yd0(x):
    return yd0 - x[3 * m]


def c_v0(x):
    return v0 - x[4 * m]


def c_vd0(x):
    return vd0 - x[5 * m]


def c_xf(x):
    return xf - x[1 * m - 1]


def c_yf(x):
    return yf - x[3 * m - 1]


con1 = {'type': 'eq', 'fun': c_x0}
con2 = {'type': 'eq', 'fun': c_xd0}
con3 = {'type': 'eq', 'fun': c_y0}
con4 = {'type': 'eq', 'fun': c_yd0}
con5 = {'type': 'eq', 'fun': c_v0}
con6 = {'type': 'eq', 'fun': c_vd0}
con7 = {'type': 'eq', 'fun': c_xf}
con8 = {'type': 'eq', 'fun': c_y0}

boundry_conditions = [con1, con2, con3, con4, con5, con6, con7, con8]

cons_total = collocation_cons + boundry_conditions


bound_pos_x = (0.0, xf)
bound_pos_y = (-10.0, 0.0)
bound_rate_x = (0.0, 10.0)
bound_rate_y = (-10.0, 0.0)
bound_rate = (0.0, 10.0)
bound_control = (-ma.pi / 2., ma.pi / 2.)
bound_time = (0.0, 2.0)
bound_list = [bound_pos_x, bound_rate_x, bound_pos_y, bound_rate_y, bound_rate, bound_rate, bound_control]

bnds = cf.bound_gen(bound_list, ns) + (bound_time,)
print(len(bnds))
for i in range(len(bnds)):
    print(bnds[i])

print("solving...")
solution = minimize(beom.objective_brach, xo, method='SLSQP',\
                    bounds=bnds, constraints=cons_total)
xp = solution.x
print(solution.nit)
print(solution.success)
print(solution.status)
print(solution.message)
print(xp)
print(len(xp))
print(beom.objective_brach(xp))