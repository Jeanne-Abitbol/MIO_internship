import matplotlib.pyplot as plt
from speeds_f import *

Tmax = 24*3
Zmax = 300
dt = 0.001
dz = 10
vd_max = 100
vr_max = 80
delta = 10
vmax = 11698
treshold = 1000

tt = np.int(Tmax / dt)
zz = np.int(Zmax / dz)
time = np.arange(0, Tmax, dt)
watercolumn = np.arange(0, Zmax, dz)

Z0 = np.zeros((tt, zz))
for i in range(zz):
    Z0[0, i] = i * dz * np.exp(-0.05 * i * dz)


P0 = np.zeros(zz)
for i in range(zz):
    P0[i] = norm(i * dz, 25, 20)

speed = v_richards_dI0dt
args = (vd_max, vr_max, delta)

s, Z = model_sweby_migrationonly(Z0, P0, dz, dt, speed, args)

plt.figure()
plt.title('speed')
plt.plot(time, s[:, 0], label='z = 0')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(200 / dz)], label='z = 200')
#plt.plot(time, s[:, np.int(500 / dz)], label='z = 500')
plt.vlines([6, 18], -5, 5, linestyles='dashed')
plt.legend()

y_max = np.max(Z)

indices = np.arange(np.int(tt * dt * 4))
plt.figure()
plt.xlim(0, Zmax)
plt.ylim(0, y_max)
Zh = []
for i in indices:
    Zh.append(Z[np.int(i / (4 * dt))])
    plt.clf()
    plt.title('t = ' + str(i / 4 % 24))
    plt.plot(watercolumn, Zh[-1])
    plt.pause(0.01)
plt.show()

from pylab import meshgrid, cm, imshow, colorbar

plt.figure()
x = np.copy(indices)
y = np.arange(0, Zmax, dz)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(Zh, cmap=cm.RdBu, aspect='auto')  # drawing the function
colorbar(im)  # adding the colobar on the right

Ztot = np.zeros(tt)
for t in range(tt):
    Ztot[t] = np.sum(Z[t])

plt.figure()
plt.title('total biomass')
plt.plot(time, Ztot)

plt.show()