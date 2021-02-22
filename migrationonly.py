import numpy as np
import matplotlib.pyplot as plt
from pylab import meshgrid, cm, imshow, colorbar

dz = 10  # depth step (m)
dt = 0.01  # time step (h)
Zmax = 500  # m
Tmax = 24 * 3  # h


def I0_richards(t, Is=1e-6, eps=1e-3):
    return Is + (1 - Is) / 2 * (
            1 + np.sin(np.pi / 12 * (t - 6)) + np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2) - np.sqrt(eps + 1))


def I_richards(t, z, gamma=0.05):
    return I0_richards(t) * np.exp(- gamma * z)


def transport_migrationonly(R0, C0, dz, dt, speed, args, schema):
    tt, zz = np.shape(C0)
    C = np.copy(C0)
    R = np.copy(R0)
    s = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        a = np.zeros((zz))
        b = np.zeros((zz))
        c = np.zeros((zz))
        for i in range(zz):
            s[t, i] = speed(t * dt, i * dz, R[i], *args)
            nu = s[t, i] * dt / dz
            # CFL condition
            print('CFL = '+str(nu))
            if np.abs(nu) >= 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None
            if schema == 'SC':
                a[i], b[i], c[i] = -nu / 2, 1, nu / 2
            elif schema == 'SLF':
                a[i], b[i], c[i] = (1 - nu) / 2, 0, (1 + nu) / 2
            elif schema == 'SDA':
                if s[t, i] >= 0:
                    a[i], b[i], c[i] = 0, 1 - nu, nu
                else:
                    a[i], b[i], c[i] = -nu, 1 + nu, 0
            elif schema == 'SLW':
                a[i], b[i], c[i] = nu * (nu - 1) / 2, 1 - nu ** 2, nu * (nu + 1) / 2

        for i in range(1, zz - 1):
            C[t + 1, i] = a[i] * C[t, i + 1] + b[i] * C[t, i] + c[i] * C[t, i - 1]
        C[t + 1, 0] = C[t + 1, 2]
        C[t + 1, -1] = C[t + 1, -3]

    return s, C


def transport2_migrationonly(R0, C0, dz, dt, speed, args):
    tt, zz = np.shape(C0)
    C = np.copy(C0)
    R = np.copy(R0)
    s = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        for i in range(zz):
            s[t, i] = speed(t * dt, i * dz, R[i], *args)
            nu = s[t, i] * dt / dz
            # CFL condition
            if np.abs(nu) >= 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

        for i in range(1, zz - 1):
            C[t + 1, i] = C[t, i] - dt / dz * ((s[t, i] > 0) * s[t, i] * (C[t, i] - C[t, i - 1]) - (s[t, i] < 0) *
                                               s[t, i] * (C[t, i + 1] - C[t, i]))
        C[t + 1, 0] = 0
        C[t + 1, -1] = 0

    return s, C


def v(t, z, R):
    th = t % 24
    if 4.5 <= th <= 6.5:
        return 72
    elif 16.5 <= th <= 18.5:
        return -72
    else:
        return 0


tt = np.int(Tmax / dt)
zz = np.int(Zmax / dz)

time = np.arange(0, Tmax, dt)
watercolumn = np.arange(0, Zmax, dz)

C0 = np.zeros((tt, zz))
for i in range(4, zz):
    C0[0, i] = (i - 4) * dz * np.exp(-0.05 * (i - 4) * dz)

R0 = np.zeros((zz))
for i in range(zz):
    R0[i] = i * dz * np.exp(-0.05 * i * dz)


light = I_richards
speed = v
args = ()
s, C = transport2_migrationonly(R0, C0, dz, dt, speed, args)

plt.figure()
plt.title('speed')
plt.plot(time, s[:, np.int(20 / dz)], label='z = 20')
plt.plot(time, s[:, np.int(350 / dz)], ':', label='z = 350')
plt.legend()

y_max = np.max(C)
y_min = np.min(C)

indices = np.arange(np.int(tt * dt * 4))
plt.figure()
Ch = []
for i in indices:
    Ch.append(C[np.int(i / (4 * dt))])
    plt.clf()
    plt.xlim(0, Zmax)
    plt.ylim(y_min, y_max)
    plt.title('t = ' + str(i / 4 % 24))
    plt.plot(watercolumn, Ch[-1])
    plt.pause(0.01)


Ctot = np.zeros(tt)
for t in range(tt):
    Ctot[t] = np.sum(C[t])

plt.figure()
plt.title('total biomass')
plt.plot(time, Ctot)
plt.legend()

plt.figure()
x = np.copy(indices)
y = np.arange(0, Zmax, dz)
X, Y = meshgrid(x, y)  # grid of point
im = imshow(Ch, cmap=cm.RdBu, aspect='auto')  # drawing the function
colorbar(im)  # adding the colobar on the right


plt.show()
