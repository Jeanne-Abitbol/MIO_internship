import matplotlib.pyplot as plt
import datetime
from pylab import meshgrid, cm, imshow, colorbar
from model_G02_f import *

currentDT = datetime.datetime.now()
start_time = 60 * currentDT.hour + currentDT.minute + currentDT.second / 60  # minutes

dz = 10  # depth step (m)
dt = 0.001  # time step (h)
Zmax = 500  # m
Tmax = 24 * 7  # h
K = 3
alpha = 0.01
beta = 1
e = 0.8
mu = 1e-2
r_max = 0.01
K_I = 0.2
delta = 10
tdI = 3
v = 250
d = 0.1

tt = np.int(Tmax / dt)
zz = np.int(Zmax / dz)

time = np.arange(0, Tmax, dt)
watercolumn = np.arange(0, Zmax, dz)

C0 = np.zeros((tt, zz))
for i in range(4, zz):
    C0[0, i] = (i - 4) * dz * np.exp(-0.05 * (i - 4) * dz)

R0 = np.zeros((tt, zz))
for i in range(zz):
    R0[0, i] = i * dz * np.exp(-0.05 * i * dz)


light = I_richards
speed = v_richards_seuil
args = (dt, tdI, v)
s, R, G, C, D = model_AN_RGCD(R0, C0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu, d)

plt.figure()
plt.title('speed')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(350 / dz)], ':', label='z = 350')
plt.legend()

y_max = np.max([np.max(R), np.max(C), np.max(D)])

indices = np.arange(np.int(tt * dt * 4))
plt.figure()
plt.xlim(0, Zmax)
plt.ylim(0, y_max)
Ch = []
for i in indices:
    Ch.append(C[np.int(i / (4 * dt))])
    plt.clf()
    plt.title('t = ' + str(i / 4 % 24))
    plt.plot(watercolumn, R[np.int(i / (4 * dt))], label='R')
    plt.plot(watercolumn, Ch[-1], label='C')
    plt.plot(watercolumn, D[np.int(i / (4 * dt))], label='D')
    plt.legend()
    plt.pause(0.01)
plt.show()

"""plt.figure()
for i in indices:
    plt.clf()
    plt.title('t = ' + str(i / 4 % 24))
    plt.plot(watercolumn, G[np.int(i / (4 * dt))], label='G')
    plt.legend()
    plt.pause(0.01)
plt.show()"""

Ctot = np.zeros(tt)
Rtot = np.zeros(tt)
Dtot = np.zeros(tt)
for t in range(tt):
    Ctot[t] = np.sum(C[t])
    Rtot[t] = np.sum(R[t])
    Dtot[t] = np.sum(D[t])

plt.figure()
plt.title('total biomass')
plt.plot(time, Ctot, label='C')
plt.plot(time, Rtot, label='R')
plt.plot(time, Dtot, label='D')
plt.legend()

Gtot = np.zeros(tt)
for t in range(tt):
    Gtot[t] = np.sum(G[t])

plt.figure()
plt.title('total gut content')
plt.plot(time, Gtot, label='G')
plt.legend()


plt.figure()
x = np.copy(indices)
y = np.arange(0, Zmax, dz)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(Ch, cmap=cm.RdBu, aspect='auto')  # drawing the function
colorbar(im)  # adding the colobar on the right

currentDT = datetime.datetime.now()
final_time = 60 * currentDT.hour + currentDT.minute + currentDT.second / 60  # minutes
running_time = final_time - start_time
print('The execution took '+str(running_time)+' minutes.')

plt.show()