import matplotlib.pyplot as plt
import datetime
from pylab import meshgrid, cm, imshow, colorbar
from model_global01_f import *

currentDT = datetime.datetime.now()
start_time = 60 * currentDT.hour + currentDT.minute + currentDT.second / 60  # minutes

dz = 10  # depth step (m)
dt = 0.001  # time step (h)
Zmax = 300  # m
Tmax = 24 * 7  # h
K = 3
alpha = 0.01
beta = 1
e = 0.8
mu = 1e-2
r_max = 0.5
K_I = 0.5
delta = 10
vr_max = 63
vd_max = 72
rho = 0.01
Em = 0.5
tdI = 0.2
td = 50
tl = 1e-4
dl = 9e-5

tt = np.int(Tmax / dt)
zz = np.int(Zmax / dz)

time = np.arange(0, Tmax, dt)
watercolumn = np.arange(0, Zmax, dz)

Z0 = np.zeros((tt, zz))
for i in range(zz):
    Z0[0, i] = i * dz * np.exp(-0.05 * i * dz)

P0 = np.zeros((tt, zz))
for i in range(zz):
    P0[0, i] = i * dz * np.exp(-0.05 * i * dz)

E0 = np.zeros((tt, zz))
for i in range(zz):
    E0[0, i] = (Z0[0, i] > 0) * Em


light = I_richards
speed = v_richards_dI0dt
args = (vd_max, vr_max, delta, tdI, td, tl, dl)
s, P, E, Z, D = model_AN_RC(P0, Z0, E0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu, rho, Em)

plt.figure()
plt.title('speed')
plt.plot(time, s[:, 0], label='z = 0')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(200 / dz)], label='z = 200')
plt.vlines([0, 4, 12, 20, 24, 28, 36, 44, 48], -5, 5, linestyles='dashed')
plt.legend()

y_max = np.max([np.max(P), np.max(Z), np.max(D)])

indices = np.arange(np.int(tt * dt * 4))
plt.figure()
plt.xlim(0, Zmax)
plt.ylim(0, y_max)
Zh = []
for i in indices:
    Zh.append(Z[np.int(i / (4 * dt))])
    plt.clf()
    plt.title('t = ' + str(i / 4 % 24))
    plt.plot(watercolumn, P[np.int(i / (4 * dt))], label='P')
    plt.plot(watercolumn, Zh[-1], label='Z')
    plt.plot(watercolumn, D[np.int(i / (4 * dt))], label='D')
    plt.legend()
    plt.pause(0.01)
plt.show()

Ztot = np.zeros(tt)
Ptot = np.zeros(tt)
Dtot = np.zeros(tt)
for t in range(tt):
    Ztot[t] = np.sum(Z[t])
    Ptot[t] = np.sum(P[t])
    Dtot[t] = np.sum(D[t])

plt.figure()
plt.title('total biomass')
plt.plot(time, Ztot, label='Z')
plt.plot(time, Ptot, label='P')
plt.plot(time, Dtot, label='D')
plt.legend()


plt.figure()
x = np.copy(indices)
y = np.arange(0, Zmax, dz)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(Zh, cmap=cm.RdBu, aspect='auto')  # drawing the function
colorbar(im)  # adding the colobar on the right

plt.show()

currentDT = datetime.datetime.now()
final_time = 60 * currentDT.hour + currentDT.minute + currentDT.second / 60  # minutes
running_time = final_time - start_time
print('The execution took '+str(running_time)+' minutes.')