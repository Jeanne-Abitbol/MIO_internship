import numpy as np
import matplotlib.pyplot as plt

"""def sr(t):
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
    print(i)"""

import matplotlib
matplotlib.use('Qt5Agg') #use Qt5 as backend, comment this line for default backend

from matplotlib import pyplot as plt
from matplotlib import animation

fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(0, 100))

N = 4
lines = [plt.plot([], [])[0] for _ in range(N)] #lines to animate

rectangles = plt.bar([0.5,1,1.5],[50,40,90],width=0.1) #rectangles to animate

patches = lines + list(rectangles) #things to animate

def init():
    #init lines
    for line in lines:
        line.set_data([], [])

    #init rectangles
    for rectangle in rectangles:
        rectangle.set_height(0)

    return patches #return everything that must be updated

def animate(i):
    #animate lines
    for j,line in enumerate(lines):
        line.set_data([0, 2], [10 * j,i])

    #animate rectangles
    for j,rectangle in enumerate(rectangles):
        rectangle.set_height(i/(j+1))

    return patches #return everything that must be updated

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True)

plt.show()
