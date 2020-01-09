import math as ma
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.gridspec as gridspec
from problem_functions import brach_eom as beom
from utilities import collocation_functions as cf

"""soln_tol of 1e-4 and ns=15 works well for xf,yf=1.0,0.3"""
soln_tol = 0.000005

ns = 15 # number of segments
m = ns + 1 # number of gridpoints
p = 3 # number of states
q = 1 # number of control variables
gamma = p + q

# Boundry Conditions / Constraints
x0 = 0.0
xf = 1.0
y0 = 0.0
yf = 0.3
v0 = 0.0

plot_initialation = False
xo = beom.ramp_initializer(x0, xf, y0, yf, v0, m, plot_initialation)

hk = xo[-1] / ns

print("xo:\n{}".format(xo))
print("len xo:\n {}".format(len(xo)))

collocation_cons = cf.constraint_list_gen(ns, p, q, cf.dfct_eval, beom.brach_eom_list)

print("len collocation cons:\n {}".format(len(collocation_cons)))

# constraints
def c_x0(x):
    return x[0 * m] - 0.0


def c_y0(x):
    return x[1 * m] - 0.0


def c_v0(x):
    return x[2 * m] - 0.0


def c_xf(x):
    return x[1 * m - 1] - xf


def c_yf(x):
    return x[2 * m - 1] - yf


con1 = {'type': 'eq', 'fun': c_x0}
con2 = {'type': 'eq', 'fun': c_y0}
con3 = {'type': 'eq', 'fun': c_v0}
con4 = {'type': 'eq', 'fun': c_xf}
con5 = {'type': 'eq', 'fun': c_yf}

boundry_conditions = [con1, con2, con3, con4, con5]

cons_total = collocation_cons + boundry_conditions

print("number of  constraints:\n{}".format(len(cons_total)))
print("all constraints:\n{}".format(cons_total))

for i in range(len(cons_total)):
    print("constraint {}\n{}".format(i+1,cons_total[i]))

bound_pos_x = (0.0, xf)
bound_pos_y = (0.0, 3.0)
bound_rate = (0.0, 6.0)
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
                    bounds=bnds, constraints=cons_total, options={"disp":True, "ftol":soln_tol})
xp = solution.x

xs = xp[0:m]
ys = xp[m:2 * m]
vs = xp[2 * m:3 * m]
us = xp[3 * m:4 * m]

print("\nx:\n {}\n".format(xs))
print("y:\n {}\n".format(ys))
print("v:\n {}\n".format(vs))
print("u:\n {}\n".format(us))
print(len(xp))
print(beom.objective_brach(xp))

optimal_states, optimal_controls = cf.state_parser(xp, m, p)
x_opt = optimal_states[0]
y_opt = optimal_states[1]
v_opt = optimal_states[2]
u_opt = optimal_controls[0]

path_opt = InterpolatedUnivariateSpline(x_opt, y_opt, k=2)
x_fine = np.linspace(0, xf, 400)


initial_states, initial_controls = cf.state_parser(xo, m, p)
x_init = initial_states[0]
y_init = initial_states[1]
v_init = initial_states[2]
u_init = initial_controls[0]

r = beom.cycloid_solver(xf, yf)[0]
theta = beom.cycloid_solver(xf, yf)[1]

xr = beom.cycloid_cordinates_xy(r, theta, 100)[0]
yr = beom.cycloid_cordinates_xy(r, theta, 100)[1]

print "\n\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
print "Travel Time on Curve Found vis Direct Collocation: ",beom.objective_brach(xp)
print "Optimal Travel Time Via Brachistocrone: ", beom.cycloid_time(xf, yf)
print "Solution via Direct Collocation is wihin {}% of Analytic Optimum.".format((1.-beom.objective_brach(xp)/beom.cycloid_time(xf, yf))*100.)
print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

fig = plt.figure()

plt.plot(x_init, -1.0*np.array(y_init), color="teal", linewidth=4.0, zorder=1, label="Initial Guess")
# plt.plot(x_opt, -1.0*np.array(y_opt), color="indianred", linewidth=4.0, zorder=3, linestyle="--", label="Optimal Solution via Direct Collocation")
plt.plot(x_fine, -1.0*np.array(path_opt(x_fine)), color="indianred", linewidth=3, linestyle=":", zorder=10, label="Optimal Solution via Direct Collocation")
plt.plot(xr, -1.0*np.array(yr), color="darkviolet", zorder=2, linewidth=4.0, label="Analyitic Optimal Solution")
plt.legend(fontsize='xx-large')
plt.grid()
plt.show()
