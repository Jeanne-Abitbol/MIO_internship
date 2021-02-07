import matplotlib.pyplot as plt
import scipy.integrate as spi
from nomigration_f import *

K=3
e=0.8
r_max = 2.32e-7
K_I = 700
alpha = 0.5
beta = 1
mu = 1e-5

Tmax = 24 * 7
time = np.arange(Tmax)

z = 0
Y0=[0.1, 0]
yap0 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# z = 20
# Y0=[0.006, 0.002]
# yap20 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# z = 250
# Y0 = [0.001, 0.01]
# yap250 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# z = 500
# Y0 = [0, 0.01]
# yap500 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

plt.figure()
plt.title('surface')
plt.plot(time, yap0[:, 0], 'b:')
plt.plot(time, yap0[:, 1], 'b')
# plt.figure()
# plt.title('z = 20')
# plt.plot(time, yap20[:, 0], 'g:')
# plt.plot(time, yap20[:, 1], 'g')
# plt.figure()
# plt.title('z = 250')
# plt.plot(time, yap250[:, 0], 'y:')
# plt.plot(time, yap250[:, 1], 'y')
# plt.figure()
# plt.title('z = 500')
# plt.plot(time, yap500[:, 0], 'r:')
# plt.plot(time, yap500[:, 1], 'r')

"""plt.figure()
plt.title('Phase plane at z = 20')
xx = np.arange(0, 0.01, 0.001)
X, Y = np.meshgrid(xx, xx)
FX = 0 * X
FY = 0 * X
n, p = np.shape(X)
z = 20
for i in range(n):
    for j in range(p):
        Z = nomigration([X[i, j], Y[i, j]], 0, z, alpha, beta, K, e, mu, r_max, K_I)
        FX[i, j] = Z[0]
        FY[i, j] = Z[1]
plt.quiver(X, Y, FX, FY)
plt.xlabel('P')
plt.ylabel('Z')"""


"""
alpha_PZ = alpha_PN = alpha_ZN = 0.5
beta_PZ = beta_PN = beta_ZN = 1
m_P = m_Z = m_N = 1e-10
mu_N = 0.01
rho_Z = rho_N = 1e-3
EZm = 0
ENm = 0

# Y = P, EZ, Z, EN, N, D
Y0 = [1, 0, 0, 0, 0, 0]

z = 0
yap0 = spi.odeint(PZNDdroop, Y0, time, args=(z, K, e, r_max, K_I, I, alpha_PZ, alpha_PN, alpha_ZN, beta_PZ, beta_PN, beta_ZN, m_P, m_Z, m_N, mu_N, rho_Z, rho_N, EZm, ENm))
print(yap0)
z = 20
yap20 = spi.odeint(PZNDdroop, Y0, time, args=(z, K, e, r_max, K_I, I, alpha_PZ, alpha_PN, alpha_ZN, beta_PZ, beta_PN, beta_ZN, m_P, m_Z, m_N, mu_N, rho_Z, rho_N, EZm, ENm))

z = 250
yap250 = spi.odeint(PZNDdroop, Y0, time, args=(z, K, e, r_max, K_I, I, alpha_PZ, alpha_PN, alpha_ZN, beta_PZ, beta_PN, beta_ZN, m_P, m_Z, m_N, mu_N, rho_Z, rho_N, EZm, ENm))

plt.figure()
plt.title('z = 0')
plt.plot(time, yap0[:, 0], label='P')
plt.plot(time, yap0[:, 2], label='Z')
plt.plot(time, yap0[:, 4], label='N')
plt.plot(time, yap0[:, 5], label='D')
plt.legend()

plt.figure()
plt.title('z = 20')
plt.plot(time, yap20[:, 0], label='P')
plt.plot(time, yap20[:, 2], label='Z')
plt.plot(time, yap20[:, 4], label='N')
plt.plot(time, yap20[:, 5], label='D')
plt.legend()

plt.figure()
plt.title('z = 250')
plt.plot(time, yap250[:,0], label='P')
plt.plot(time, yap250[:,2], label='Z')
plt.plot(time, yap250[:,4], label='N')
plt.plot(time, yap250[:,5], label='D')
plt.legend()"""

plt.show()

