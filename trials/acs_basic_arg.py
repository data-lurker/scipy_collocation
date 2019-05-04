import numpy as np
from scipy.optimize import minimize

ns = 5
M = ns + 1

I = 1000.0


def xd1(x2):
    return x2


def xd2(u):
    return u / I


# ykpl := "y k plus 1"
def yc(yk, ykp1, fk, fkp1, hk):
    return .5 * (yk + ykp1) + (hk / 8.0) * (fk - fkp1)


def defect(yk, ykp1, fk, fkp1, fc, hk):
    return ykp1 - yk - (hk / 6.0) * (fk + 4.0 * fc + fkp1)


t1 = 0.0
tf = 5.0
hk = (tf - t1) / ns
t = np.linspace(t1, tf, 100)

# nonlinear programme variables, initial guesse
y1 = np.linspace(10.0, 0.0, M) + np.random.uniform(-3.0, 3.0, M)
y1[0] = 10.0
y2 = np.random.uniform(-3.0, 3.0, M)
u = np.random.uniform(-5.0, 5.0, M)
xo = np.concatenate((y1, y2, u), axis=0)

# constraints
xo[0] = 10.0
xo[M - 1] = 0.0
xo[M] = 0.0
xo[2 * M - 1] = 0.0


def c1(x):
    return x[0] - 10.0


def c2(x):
    return x[M - 1]


def c3(x):
    return x[M]


def c4(x):
    return x[2 * M - 1]

# ======================================================
# segment 1
def da1(x, n):
    yk = x[n - 1]
    ykp1 = x[n]
    fk = xd1(x[M + n - 1])
    fkp1 = xd1(x[M + n])
    ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd1(ykc)
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def db1(x, n):
    yk = x[M + n - 1]
    ykp1 = x[M + n]
    fk = xd2(x[2 * M + n - 1])
    fkp1 = xd2(x[2 * M + n])
    # ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd2(xd2((x[2 * M + n - 1] + x[2 * M + n]) / 2.0))
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d
# ======================================================
# segment 2
def da2(x, n):
    yk = x[n - 1]
    ykp1 = x[n]
    fk = xd1(x[M + n - 1])
    fkp1 = xd1(x[M + n])
    ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd1(ykc)
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def db2(x, n):
    yk = x[M + n - 1]
    ykp1 = x[M + n]
    fk = xd2(x[2 * M + n - 1])
    fkp1 = xd2(x[2 * M + n])
    # ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd2(xd2((x[2 * M + n - 1] + x[2 * M + n]) / 2.0))
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d
# ======================================================
# segment 3
def da3(x, n):
    yk = x[n - 1]
    ykp1 = x[n]
    fk = xd1(x[M + n - 1])
    fkp1 = xd1(x[M + n])
    ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd1(ykc)
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def db3(x, n):
    yk = x[M + n - 1]
    ykp1 = x[M + n]
    fk = xd2(x[2 * M + n - 1])
    fkp1 = xd2(x[2 * M + n])
    # ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd2(xd2((x[2 * M + n - 1] + x[2 * M + n]) / 2.0))
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d
# ======================================================
# segment 4
def da4(x, n):
    yk = x[n - 1]
    ykp1 = x[n]
    fk = xd1(x[M + n - 1])
    fkp1 = xd1(x[M + n])
    ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd1(ykc)
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def db4(x, n):
    yk = x[M + n - 1]
    ykp1 = x[M + n]
    fk = xd2(x[2 * M + n - 1])
    fkp1 = xd2(x[2 * M + n])
    # ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd2(xd2((x[2 * M + n - 1] + x[2 * M + n]) / 2.0))
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d
# ======================================================
# segment 5
def da5(x, n):
    yk = x[n - 1]
    ykp1 = x[n]
    fk = xd1(x[M + n - 1])
    fkp1 = xd1(x[M + n])
    ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd1(ykc)
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def db5(x, n):
    yk = x[M + n - 1]
    ykp1 = x[M + n]
    fk = xd2(x[2 * M + n - 1])
    fkp1 = xd2(x[2 * M + n])
    # ykc = yc(yk, ykp1, fk, fkp1, hk)
    fc = xd2(xd2((x[2 * M + n - 1] + x[2 * M + n]) / 2.0))
    d = defect(yk, ykp1, fk, fkp1, fc, hk)
    return d


def objective(x):
    tot = 0.
    q11 = 4.0
    R = 1.0
    for i in range(M):
        tot += q11 * x[i]**2
        tot += R * x[i + 2 * M]**2
    return tot


bound_pos = (-50.0, 50.0)
bound_rate = (-1000.0, 1000.0)
bound_control = (-200.0,200.0)
bnds = (bound_pos, bound_pos, bound_pos, bound_pos, bound_pos,bound_pos,
        bound_rate, bound_rate, bound_rate, bound_rate, bound_rate, bound_rate,
        bound_control, bound_control, bound_control, bound_control, bound_control, bound_control)

con1 = {'type': 'eq', 'fun': c1}
con2 = {'type': 'eq', 'fun': c2}
con3 = {'type': 'eq', 'fun': c3}
con4 = {'type': 'eq', 'fun': c4}

con_da1 = {'type': 'eq', 'fun': da1, 'args': (1,)}
con_db1 = {'type': 'eq', 'fun': db1, 'args': (1,)}

con_da2 = {'type': 'eq', 'fun': da2, 'args': (2,)}
con_db2 = {'type': 'eq', 'fun': db2, 'args': (2,)}

con_da3 = {'type': 'eq', 'fun': da3, 'args': (3,)}
con_db3 = {'type': 'eq', 'fun': db3, 'args': (3,)}

con_da4 = {'type': 'eq', 'fun': da4, 'args': (4,)}
con_db4 = {'type': 'eq', 'fun': db4, 'args': (4,)}

con_da5 = {'type': 'eq', 'fun': da5, 'args': (5,)}
con_db5 = {'type': 'eq', 'fun': db5, 'args': (5,)}

cons = ([con1, con2, con3,con4,
        con_da1, con_db1,
        con_da2, con_db2,
        con_da3, con_db3,
        con_da4, con_db4,
        con_da5, con_db5,
       ])

print("ready...")
print("origional_objective: {}".format(objective(xo)))
solution = minimize(objective, xo, method='SLSQP',\
                    bounds=bnds, constraints=cons)
x = solution.x

print(x)
print(x[0:M-1])
print(x[M:M-1])
print(x[(2*M):(2*M-1)])
print('Final Objective: ' + str(objective(x)))
print("done.")

# TODO create overlay of optimal points vs. recovered points
# TODO integarte optimal points, compare to third order spline fit

