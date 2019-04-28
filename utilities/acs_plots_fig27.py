import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from integrators import euler 

path_x1 = r"/home/rbyron/Documents/research/python_projekts/optimization/data/x1.csv"
path_x2 = r"/home/rbyron/Documents/research/python_projekts/optimization/data/x2.csv"
path_u = r"/home/rbyron/Documents/research/python_projekts/optimization/data/u.csv"

x1 = np.loadtxt(open(path_x1, "r"), delimiter=",")
x2 = np.loadtxt(open(path_x2, "r"), delimiter=",")
u = np.loadtxt(open(path_u, "r"), delimiter=",")

x = np.linspace(0.0, 15.0, 100)

spl_x1 = UnivariateSpline([x1[i][0] for i in range(len(x1))], [x1[i][1] for i in range(len(x1))], k=3, ext=0, s=0)
spl_x2 = UnivariateSpline([x2[i][0] for i in range(len(x2))], [x2[i][1] for i in range(len(x2))], k=3, ext=0, s=0)
spl_u = UnivariateSpline([u[i][0] for i in range(len(u))], [u[i][1] for i in range(len(u))], k=1, ext=0, s=0)

x1i = 10.0
x2i = 0.0

x1i_ar = []
x2i_ar = []

t0 = 0.0
tf = 12.0
h = 0.001
t = np.arange(t0, tf, h)

I = 1.0

for i in t:
    x1i = x1i + h * x2i
    x2i = x2i + h * spl_u(i) / I
    x1i_ar.append(x1i)
    x2i_ar.append(x2i)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot([x1[i][0] for i in range(len(x1))], [x1[i][1] for i in range(len(x1))], color='purple')
ax.plot([x2[i][0] for i in range(len(x2))], [x2[i][1] for i in range(len(x2))], color='aquamarine')

ax.plot(x, spl_x1(x), color='indianred')
ax.plot(x, spl_x2(x), color='goldenrod')

ax.plot(t, x1i_ar, color="cyan", linestyle="--")
ax.plot(t, x2i_ar, color="black", linestyle="--")

plt.grid()
plt.show()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot([u[i][0] for i in range(len(u))], [u[i][1] for i in range(len(u))], color='purple')

ax.plot(x, spl_u(x), color='indianred')

plt.grid()
plt.show()
