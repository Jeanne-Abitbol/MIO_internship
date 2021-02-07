import numpy as np
import matplotlib.pyplot as plt

def sr(t):
    return np.sin(2*np.pi*t/16+np.pi/2)


Tmax = 20
dt = 0.5
time = np.arange(0, 72, dt)
plt.figure()
plt.plot(time, sr(time))
plt.vlines([4, 12, 20], -1, 1, linestyles='dashed')
plt.hlines(0, 0, 24, linestyles='dashed')


def I0_poggiale(t, alpha):
    th = t % 24
    return 1500 * np.exp((-alpha*(th-12)**2)/((th-12)**2-64)**2)

plt.figure()
plt.title('Poggiale light')
plt.plot(time, I0_poggiale(time, 1))


def norm(x, mu, sigma):
    return np.exp(-0.5*((x-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi))

time = np.arange(600)
plt.figure()
plt.plot(time, norm(time, 30, 20))

plt.show()

for i in(2,3):
    print(i)