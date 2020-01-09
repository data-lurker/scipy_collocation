import matplotlib.pyplot as plt
from brach_eom import ramp_initializer, slide_down_ramp, cycloid_solver, cycloid_cordinates_xy

x0 = 0.0
xf = 1.0
y0 = 0.0
yf = 0.3
v0 = 0.0
hk = 0.00001
m = 5

# state_history = ramp_initializer(x0, xf, y0, yf, v0, m, False)

# state_history = slide_down_ramp(x0, xf, y0, yf, v0, hk)

# print state_history

r = cycloid_solver(xf, yf)[0]
theta = cycloid_solver(xf, yf)[1]

xr = cycloid_cordinates_xy(r, theta, 100)[0]
yr = cycloid_cordinates_xy(r, theta, 100)[1]
plt. plot(xr, yr)
plt.show()


# print("x position history: {}".format(state_history[0]))
# print("y position history: {}".format(state_history[1]))
# print("path velocity history: {}".format(state_history[2]))
# print("Final Time: {}".format(state_history[3]))
