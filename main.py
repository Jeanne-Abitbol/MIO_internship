import matplotlib.pyplot as plt
from model_functions import *
import datetime


currentDT = datetime.datetime.now()
start_time = 60 * currentDT.hour + currentDT.minute + currentDT.second / 60  # minutes

dz = 4  # depth step (m)
dt = 0.001  # time step (h)
Zmax = 600  # m
Tmax = 24 * 3  # h
K = 3
alpha = 0.5
beta = 1
e = 0.8
mu = 0
r_max = 1e-6
K_I = 700
delta = 10
vr_max = 72
vd_max = 72

tt = np.int(Tmax / dt)
zz = np.int(Zmax / dz)

time = np.arange(0, Tmax, dt)
watercolumn = np.arange(0, Zmax, dz)

Z0 = np.zeros((tt, zz))
for i in range(zz):
    Z0[0, i] = norm(i * dz, 50, 10)

# P0 = np.zeros((zz))
# for i in range(zz):
#    P0[i] = norm(i * dz, 25, 20)

P0 = np.zeros((tt, zz))
for i in range(zz):
    P0[0, i] = norm(i * dz, 25, 20)

"""s, Z = model_AN_migrationonly(Z0, P0, dz, dt, v_madani2, (vd_max, vr_max, delta))
# print(Z)


plt.figure()
plt.title('speed')
plt.plot(time, s[:, 0], label='z = 0')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(150 / dz)], label='z = 150')
plt.plot(time, s[:, np.int(500 / dz)], label='z = 500')
plt.vlines([0, 4, 12, 20, 24, 28, 36, 44, 48], -100, 100, linestyles='dashed')
plt.legend()

plt.figure()
for i in range(np.int(tt * dt * 4)):
    plt.clf()
    plt.title('t = ' + str(i / 4))
    plt.plot(watercolumn, Z[np.int(i / (4 * dt))])
    plt.pause(0.01)
plt.show()

Ztot = np.zeros(tt)
for t in range(tt):
    Ztot[t] = np.sum(Z[t])
plt.figure()
plt.title('total zooplankton biomass')
plt.plot(time, Ztot)
# print(Ztot)"""

"""maxZ = np.zeros((tt))
for t in range(tt):
    maxZ[t] = np.argmax(Z[t, :])

# days = np.arange(7)
# days *= 24

plt.figure()
plt.title('maxZ')
plt.plot(time, maxZ)
# plt.vlines(days, 0, np.max(maxZ), linestyles='dashed')
plt.show()"""

# plt.figure()
# plt.plot(watercolumn, P0[0], label='P0')
# plt.plot(watercolumn, Z0[0], label='Z0')
# plt.legend()


"""s, P, Z = model_AN(Z0, P0, dz, dt, I, v_madani2, (vd_max, vr_max, delta), K_I, r_max, alpha, beta, K, e, mu)


plt.figure()
plt.title('speed')
plt.plot(time, s[:, 0], label='z = 0')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(150 / dz)], label='z = 150')
#plt.plot(time, s[:, np.int(500 / dz)], label='z = 500')
plt.vlines([0, 4, 12, 20, 24, 28, 36, 44, 48], -5, 5, linestyles='dashed')
plt.legend()

plt.figure()
for i in range(np.int(tt * dt * 4)):
    plt.clf()
    plt.title('t = ' + str(i / 4 % 24))
    plt.plot(watercolumn, P[np.int(i / (4 * dt))], label='P')
    plt.plot(watercolumn, Z[np.int(i / (4 * dt))], label='Z')
    plt.legend()
    plt.pause(0.01)
plt.show()

Ztot = np.zeros(tt)
Ptot = np.zeros(tt)
for t in range(tt):
    Ztot[t] = np.sum(Z[t])
    Ptot[t] = np.sum(P[t])
plt.figure()
plt.title('total biomass')
plt.plot(time, Ztot, label='Z')
plt.plot(time, Ptot, label='P')
plt.legend()"""

"""from pylab import meshgrid, cm, imshow, colorbar

plt.figure()
x = np.arange(0, Tmax)
y = np.arange(0, Zmax, dz)
X, Y = meshgrid(x, y)  # grid of point
hours = np.arange(Tmax)
Zh = Z[hours]
print(np.shape(Zh))
im = imshow(Zh, cmap=cm.RdBu)  # drawing the function
colorbar(im)  # adding the colobar on the right"""

rho = 1
Em = 1
E0 = np.zeros((tt, zz))
for i in range(zz):
    E0[0, i] = (Z0[0, i] > 0) * Em
func = model_AN_RC
args = (Z0, P0, E0, dz, dt, I, v_madani2_relative, (vd_max, vr_max, delta), K_I, r_max, alpha, beta, K, e, mu, rho, Em)
s, P, Z, D = model(func, args)

plt.figure()
plt.title('speed')
plt.plot(time, s[:, 0], label='z = 0')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(150 / dz)], label='z = 150')
# plt.plot(time, s[:, np.int(500 / dz)], label='z = 500')
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

from pylab import meshgrid, cm, imshow, colorbar

plt.figure()
x = np.copy(indices)
y = np.arange(0, Zmax, dz)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(Zh, cmap=cm.RdBu)  # drawing the function
colorbar(im)  # adding the colobar on the right

plt.show()

currentDT = datetime.datetime.now()
final_time = 60 * currentDT.hour + currentDT.minute + currentDT.second / 60  # minutes
running_time = final_time - start_time
print('The execution took '+str(running_time)+' minutes.quit')