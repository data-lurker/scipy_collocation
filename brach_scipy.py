import math as ma
import numpy as np
from scipy.optimize import minimize
from problem_functions import brach_eom as beom
from utilities import collocation_functions as cf

ns = 8
m = ns + 1
p = 3
q = 1
gamma = p + q

# initialize
t0 = 0.0
tf = 1.5
hk = (tf - t0) / ns

# Boundry Conditions / Constraints
# Initial
x0 = 0.0
y0 = 0.0
v0 = 0.0
# t0 = 0.0

# Final
xf = 1.0
yf = .3


w = np.linspace(x0, xf, m)
w[-1] = xf
y = np.linspace(y0, yf, m)
y[-1] = yf
v = np.linspace(v0, 5., m)
u = np.ones([m]) * np.tan((yf - y0) / (xf - x0))
# np.random.uniform(-.05, .05, m)

xo = np.concatenate((w, y, v, u, [tf]), axis=0).tolist()

print(xo)
print(len(xo))

collocation_cons = cf.constraint_list_gen(ns, p, q, cf.dfct_eval, beom.brach_eom_list)

print(len(collocation_cons))
# for i in range(len(collocation_cons)):
#     print(collocation_cons[i])


# constraints
def c_x0(x):
    return x[0 * m] - 0.0


def c_y0(x):
    return x[1 * m] - 0.0


def c_v0(x):
    return x[2 * m] - 0.0


def c_xf(x):
    return x[1 * m - 1] - 1.0


def c_yf(x):
    return x[2 * m - 1] - .3


con1 = {'type': 'eq', 'fun': c_x0}
con2 = {'type': 'eq', 'fun': c_y0}
con3 = {'type': 'eq', 'fun': c_v0}
con4 = {'type': 'eq', 'fun': c_xf}
con5 = {'type': 'eq', 'fun': c_yf}

boundry_conditions = [con1, con2, con3, con4, con5]

cons_total = collocation_cons + boundry_conditions

for i in range(len(cons_total)):
    print(cons_total[i])


bound_pos_x = (0.0, 10.0)
bound_pos_y = (0.0, 10.0)
bound_rate = (0.0, 11.0)
bound_control = (0.0, ma.pi)
bound_time = (0.05, 2.0)
bound_list = [bound_pos_x, bound_pos_y, bound_rate, bound_control]

bnds = cf.bound_gen(bound_list, ns) + (bound_time,)
print(len(bnds))
for i in range(len(bnds)):
    print(bnds[i])


for i in range(len(collocation_cons)):
    arg = collocation_cons[i]['args']
    a = (xo,) + arg
    print("defect {}: {}".format(i, cf.dfct_eval(a[0], a[1], a[2], a[3], a[4], a[5], a[6])))


print("solving...")
solution = minimize(beom.objective_brach, xo, method='SLSQP',\
                    bounds=bnds, constraints=cons_total)
xp = solution.x
print(solution.nit)
print(solution.success)
print(solution.status)
print(solution.message)

xs = xp[0:m]
ys = xp[m:2 * m]
vs = xp[2 * m:3 * m]
us = xp[3 * m:4 * m]

print("\nx: {}\n".format(xs))
print("y: {}\n".format(ys))
print("v: {}\n".format(vs))
print("u: {}\n".format(us))
print(len(xp))
print(beom.objective_brach(xp))
