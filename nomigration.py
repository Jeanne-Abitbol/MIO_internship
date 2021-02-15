import matplotlib.pyplot as plt
import scipy.integrate as spi
from nomigration_f import *

K = 3
e = 0.8
r_max = 0.5
rho = 0.01
K_I = 700
alpha = 0.01
beta = 1
#m = 1e-2
mu = 1e-3
Em = 1


Tmax = 24 * 150
time = np.arange(Tmax)

z = 30

"""
Y0 = 1
yap = spi.odeint(phyto, Y0, time, args=(z, I, K_I, r_max, K))

Rt = np.zeros((len(time)))
for t in range(len(time)):
    Rt[t] = R(t, z, r_max, K_I, I)

plt.figure()
plt.title('Phytoplankton light dependent growth rate')
plt.plot(time, Rt)

plt.figure()
plt.title('Phytoplankton light dependent growth')
plt.plot(time, yap)"""

# z = 0
# Y0=[0.1, 0]
# yap0 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# z = 20
# Y0=[0.006, 0.002]
# yap20 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# z = 250
# Y0 = [0.001, 0.01]
# yap250 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# z = 500
# Y0 = [0, 0.01]
# yap500 = spi.odeint(PZ_nomigration, Y0, time, args=(z, alpha, beta, K, e, mu, r_max, K_I))

# plt.figure()
# plt.title('surface')
# plt.plot(time, yap0[:, 0], 'b:')
# plt.plot(time, yap0[:, 1], 'b')
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


"""Y0 = [1, 0.5*Em, 0.5]
RHO = [0, 0.01, 0.1, 1, 10]
yap = []
for rho in RHO:
    yap.append(spi.odeint(RC_droop, Y0, time, args=(z, I, K_I, r_max, alpha, beta, K, e, mu, rho, Em)))

for i in range (len(RHO)):
    plt.figure()
    plt.title('Droop with no migration - rho = '+str(RHO[i]))
    plt.plot(time, yap[i][:, 0], label='R')
    plt.plot(time, yap[i][:, 1], label='E')
    plt.plot(time, yap[i][:, 2], label='C')
    plt.legend()
    plt.savefig('droop_nomigration_rho'+str(i)+'.png', format='png')"""


"""yap = spi.odeint(RC_droop, Y0, time, args=(z, I, K_I, r_max, alpha, beta, K, e, mu, rho, Em))
plt.figure()
plt.title('Droop with no migration')
plt.plot(time, yap[:, 0], label='R')
plt.plot(time, yap[:, 1], label='E')
plt.plot(time, yap[:, 2], label='C')
plt.legend()"""

"""RHO = [0, 0.01, 0.1, 1, 10]
ALPHA = [0, 0.01, 0.1, 1, 10]
yap = np.zeros((len(ALPHA), len(RHO), len(time), 3))
for i in range(len(ALPHA)):
    alpha = ALPHA[i]
    for j in range(len(RHO)):
        rho = RHO[j]
        yap[i, j] = spi.odeint(RC_droop, Y0, time, args=(z, I, K_I, r_max, alpha, beta, K, e, mu, rho, Em))

for i in range(len(ALPHA)):
    for j in range(len(RHO)):
        plt.figure()
        plt.title('Droop with no migration - alpha = '+str(ALPHA[i])+' - rho = '+str(RHO[j]))
        plt.plot(time, yap[i, j, :, 0], label='R')
        plt.plot(time, yap[i, j, :, 1], label='E')
        plt.plot(time, yap[i, j, :, 2], label='C')
        plt.legend()
        plt.savefig('droop_nomigration_alpha'+str(i)+'rho'+str(j)+'.png', format='png')


Y0 = [1, 0.5]
yap = []
for i in range(len(ALPHA)):
    alpha = ALPHA[i]
    yap.append(spi.odeint(PZ_nomigration, Y0, time, args=(z, I, K_I, r_max, alpha, beta, K, e, mu)))
for i in range(len(ALPHA)):
    plt.figure()
    plt.title('No migration - alpha = ' + str(ALPHA[i]))
    plt.plot(time, yap[i][:, 0], label='R')
    plt.plot(time, yap[i][:, 1], label='C')
    plt.legend()
    plt.savefig('nomigration_alpha' + str(i) + '.png', format='png')"""


Y0 = [1, 0.5]
PZ = spi.odeint(PZ_nomigration, Y0, time, args=(z, I_richards, K_I, r_max, alpha, beta, K, e, mu))
Y0 = [1, 0.5*Em, 0.5]
PEZ = spi.odeint(RC_droop, Y0, time, args=(z, I_richards, K_I, r_max, alpha, beta, K, e, mu, rho, Em))
Y0 = [1, 0.5, 0]
PZD = spi.odeint(PZ_nomigration_withD, Y0, time, args=(z, I_richards, K_I, r_max, alpha, beta, K, e, mu))
Y0 = [1, 0.5*Em, 0.5, 0]
PEZD = spi.odeint(RC_droop_withD, Y0, time, args=(z, I_richards, K_I, r_max, alpha, beta, K, e, mu, rho, Em))

plt.figure()
plt.title('PZ')
plt.plot(time, PZ[:, 0], 'b', label='P')
plt.plot(time, PZ[:, 1], 'g', label='Z')
plt.legend()

plt.figure()
plt.title('PEZ')
plt.plot(time, PEZ[:, 0], 'b', label='P')
plt.plot(time, PEZ[:, 1], 'y', label='E')
plt.plot(time, PEZ[:, 2], 'g', label='Z')
plt.legend()

plt.figure()
plt.title('PZD')
plt.plot(time, PZD[:, 0], 'b', label='P')
plt.plot(time, PZD[:, 1], 'g', label='Z')
plt.plot(time, PZD[:, 2], 'r', label='D')
plt.legend()

plt.figure()
plt.title('PEZD')
plt.plot(time, PEZD[:, 0], 'b', label='P')
plt.plot(time, PEZD[:, 1], 'y', label='E')
plt.plot(time, PEZD[:, 2], 'g', label='Z')
plt.plot(time, PEZD[:, 3], 'r', label='D')
plt.legend()

plt.show()