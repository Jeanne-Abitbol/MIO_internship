# Modèle de migrations verticales

import numpy as np
import pylab as plt


def Profil0(alpha, zmax, dz):
    nz = int(zmax / dz)
    N = [];
    Prof = [];
    for i in range(nz):
        z = i * dz
        Prof = Prof + [z]
        N = N + [z * np.exp(-alpha * z)]
    return N, Prof


v = 1;
tmax = 2 * 24;
zmax = 150;
dt = 0.1;
dz = 0.2;
alpha = 0.05;

nt = int(tmax / dt);
nz = int(zmax / dz);

N0, Prof = Profil0(alpha, zmax, dz)
fig, graph0 = plt.subplots()
plt.plot(Prof, N0)
N1 = N0;
N = [N0]
F = [0. for i in range(nz)];
nu = [0. for i in range(nz)];
r = [0. for i in range(nz)]
for k in range(1, nt):
    N2 = [0 for i in range(nz)]
    F[0] = v * N1[0];
    for i in range(1, nz - 1):
        F[i] = v * N1[i];
        F[i + 1] = v * N1[i + 1]
        if N1[i] == N1[i - 1]:
            nu[i - 1] = 1
        else:
            nu[i - 1] = 1 - dt / dz * (F[i] - F[i - 1]) / (N1[i] - N1[i - 1])
        if N1[i + 1] == N1[i]:
            nu[i] = 1
        else:
            nu[i] = 1 - dt / dz * (F[i + 1] - F[i]) / (N1[i + 1] - N1[i])
        if F[i + 1] == F[i]:
            r[i] = 1
        else:
            r[i] = nu[i - 1] * (F[i] - F[i - 1]) / (nu[i] * (F[i + 1] - F[i]))
        if r[i] < 0:
            r[i] = 0
        else:
            if r[i] > 1:
                r[i] = 1
        if i > 1:
            N2[i] = N1[i] - dt / dz * (F[i] - F[i - 1]) + dt / (2 * dz) * (
                        r[i] * nu[i] * (F[i + 1] - F[i]) - r[i - 1] * nu[i - 1] * (F[i] - F[i - 1]))
        else:
            N2[1] = N1[1] - dt / dz * (F[1] - F[0]) + dt / (2 * dz) * (r[1] * nu[1] * (F[2] - F[1]))
    N2[nz - 1] = N2[nz - 3]
    N2[0] = 0
    N1 = N2
    N = N + [N1]
    plt.plot(Prof, N1)

fig, graph1 = plt.subplots()
plt.plot(Prof, N0, label='N0')
plt.legend()
plt.show()
