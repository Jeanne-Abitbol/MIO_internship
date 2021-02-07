import matplotlib.pyplot as plt
from model_functions import *

dz = 10  # depth step (m)
dt = 0.0001  # time step (h)
Zmax = 300  # m
Tmax = 72 # h
K = 3
alpha = 0.8
beta = 1
e = 0.8
mu = 0.001
r_max = 2.32e-7
K_I = 700
delta = 10
vr_max = 1
vd_max = 72
wm = v1 = 72

tt = np.int(Tmax / dt)
zz = np.int(Zmax / dz)

time = np.arange(0, Tmax, dt)
watercolumn = np.arange(0, Zmax, dz)

Z0 = np.zeros((tt, zz))
for i in range(zz):
    Z0[0, i] = norm(i * dz, 30, 20)

# P0 = np.zeros((zz))
# for i in range(zz):
#     P0[i] = norm(i * dz, 25, 20)

P0 = np.zeros((tt, zz))
for i in range(zz):
    P0[0, i] = norm(i * dz, 25, 20)

V0 = np.zeros(tt)
V20 = np.zeros(tt)
V100 = np.zeros(tt)
for t in range(tt):
    V0[t] = v_madani(time[t], 0, P0[0, 0], v1, vr_max, delta)
    V20[t] = v_madani(time[t], 20, P0[0, np.int(20/dz)], v1, vr_max, delta)
    V100[t] = v_madani(time[t], 100, P0[0, np.int(100/dz)], v1, vr_max, delta)

plt.figure()
plt.title('Madani migration speed')
plt.plot(time, V0, label='z = 0')
plt.plot(time, V20, label='z = 20')
plt.plot(time, V100, label='z = 100')
plt.legend()

VV0 = np.zeros(tt)
VV20 = np.zeros(tt)
VV100 = np.zeros(tt)
for t in range(tt):
    VV0[t] = v_madani2(time[t], 0, P0[0, 0], vd_max, vr_max, delta)
    VV20[t] = v_madani2(time[t], 20, P0[0, np.int(20/dz)], vd_max, vr_max, delta)
    VV100[t] = v_madani2(time[t], 100, P0[0, np.int(100/dz)], vd_max, vr_max, delta)

plt.figure()
plt.title('Madani2 migration speed')
plt.plot(time, VV0, label='z = 0')
plt.plot(time, VV20, label='z = 20')
plt.plot(time, VV100, label='z = 100')
plt.legend()


plt.figure()
plt.title('Richards derivative')
plt.plot(time, dIdt_richards(time, 0), label='z = 0')
plt.plot(time, dIdt_richards(time, 20), label='z = 20')
plt.plot(time, dIdt_richards(time, 200), label='z = 200')
#plt.plot(time, dIdt_richards(time, 500), label='z = 500')
plt.legend()

vt0 = np.zeros(tt)
vt20 = np.zeros(tt)
vt200 = np.zeros(tt)
vt500 = np.zeros(tt)

for t in range(tt):
    vt0[t] = v_richards(t*dt, 0, P0[0, 0], v1, delta)
    vt20[t] = v_richards(t*dt, 20, P0[0, np.int(20/dz)], v1, delta)
    vt200[t] = v_richards(t*dt, 200, P0[0, np.int(200/dz)], v1, delta)
#    vt500[t] = v_richards(t*dt, 500, P0[0, np.int(500/dz)], v1, delta)

plt.figure()
plt.title('Richards light')
plt.plot(time, I_richards(time, 0), label='z = 0')
plt.plot(time, I_richards(time, 20), label='z = 20')
plt.plot(time, I_richards(time, 200), label='z = 200')
#plt.plot(time, I_richards(time, 500), label='z = 500')
plt.legend()

plt.figure()
plt.title('Richards speed')
plt.plot(time, vt0, label='z=0')
plt.plot(time, vt20, label='z=20')
plt.plot(time, vt200, label='z=200')
#plt.plot(time, vt500, label='z=500')
plt.legend()

plt.figure()
plt.title('Diel cycle of surface light intensity (Richards et al. 1996)')
plt.plot(time, I0_richards(time))
plt.vlines([6, 18], 0, 1, linestyles='dashed')

plt.figure()
plt.title('Diel cycle of surface light intensity (Kamykowski et al. 1994)')
plt.plot(time, I0_kamykowski(time))

It = np.zeros(tt)
for j in range(tt):
    It[j] = I0(time[j])

plt.figure()
plt.title('Diel cycle of surface light intensity (Andersen and Nival 1991)')
plt.plot(time, It)

plt.figure()
plt.title('Richards absolute light and rate of change of light at z = 20')
plt.plot(time, I_richards(time, 20), label='I')
plt.plot(time, dIdt_richards(time, 20), label='dI/dt')
plt.legend()

plt.figure()
plt.title('phytoplankton growth rate')
R0 = np.zeros(tt)
R25 = np.zeros(tt)
R250 = np.zeros(tt)
R500 = np.zeros(tt)
for t in range(tt):
    R0[t] = R(t*dt, 0, r_max, K_I, I)
    R25[t] = R(t*dt, 25, r_max, K_I, I)
    R250[t] = R(t*dt, 250, r_max, K_I, I)
#    R500[t] = R(t*dt, 500, r_max, K_I, I)
plt.plot(time, R0, label='z=0')
plt.plot(time, R25, label='z=25')
plt.plot(time, R250, label='z=250')
#plt.plot(time, R500, label='z=500')
plt.legend()

plt.show()

# how to do a color map
"""from numpy import exp, arange
from pylab import meshgrid, cm, imshow, contour, clabel, colorbar, axis, title, show


# the function that I'm going to plot
def z_func(x, y):
    return (1 - (x ** 2 + y ** 3)) * exp(-(x ** 2 + y ** 2) / 2)


x = arange(-3.0, 3.0, 0.1)
y = arange(-3.0, 3.0, 0.1)
X, Y = meshgrid(x, y)  # grid of point
Z = z_func(X, Y)  # evaluation of the function on the grid

im = imshow(Z, cmap=cm.RdBu)  # drawing the function
# adding the Contour lines with labels
cset = contour(Z, arange(-1, 1.5, 0.2), linewidths=2, cmap=cm.Set2)
clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
colorbar(im)  # adding the colobar on the right
# latex fashion title
title('$z=(1-x^2+y^3) e^{-(x^2+y^2)/2}$')
show()

plt.show()"""



"""def f(a, b, c):
    return a + b + c


def g(a, b):
    return a * b


def h(func, args):
    returnfunc(1, *args)


print(h(f, (2, 3)))"""

