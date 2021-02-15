from model_functions import *


def vd_madani(t, z, v1, gamma=0.05):
    return v1 * np.exp(-gamma * z) * np.pi / 8 * np.cos(2 * np.pi * t / 16 - np.pi)  # vd = v1*dI/dt


def vr_madani(t, P, v_rmax, delta):
    return - v_rmax / (1 + delta * P) * np.sin(2 * np.pi * t / 16 + np.pi / 2)


def v_madani(t, z, P, vd_max, vr_max, delta):  # c'est nimp, période de 48h
    th = t % 24
    if 4 <= th <= 20:  # go down by day
        return vd_madani(t, z, vd_max)
    else:  # go up by night
        return vr_madani(t, P, vr_max, delta)


def v_madani2(t, z, P, vd_max, vr_max, delta, gamma=0.05):
    th = t % 24
    if th <= 12:
        if I(t, z) == 0:
            return 0
        else:
            return vd_max * np.exp(- gamma * z) * th / 8 * np.cos(np.pi * th / 8 - np.pi)
    else:
        if I(t, z) == 0:
            return 0
        else:
            return -vr_max / (1 + delta * P) * np.sin(2 * np.pi * th / 16 + np.pi / 2)


def v_madani2_relative(t, z, P, vd_max, vr_max, delta):
    return v_madani2(t, z, P, vd_max, vr_max, delta) / (1 + I(t, z) / 750)


def v_richards(t, z, P, delta):
    if dIdt_richards(t, z) <= 0:
        return dIdt_richards(t, z) / (1 + delta * P)
    else:
        return dIdt_richards(t, z)


def v_richards_relative(t, z, P, vd_max, vr_max, delta, dt):
    if dIdt_richards(t, z) >= 0:
        return vd_max / dt * (I_richards(t + dt, z) / I_richards(t, z) - 1)
    else:
        return vr_max / dt * (I_richards(t + dt, z) / I_richards(t, z) - 1) / (1 + delta * P)


def v_richards_relative_brut(t, z, P, vd_max, vr_max, vmax, delta, treshold):
    if np.abs(dIdt_richards(t, z)/I_richards(t, z)) > treshold:
        if dIdt_richards(t, z) > 0:
            return dIdt_richards(t, z) / I_richards(t, z) *(vd_max/vmax)
        else:
            return dIdt_richards(t, z) / I_richards(t, z) * (vr_max / vmax) / (1 + delta*P)
    else:
        return 0

def v_richards_dIdt_noz(t, z, P, vd_max, vr_max, delta, eps = 1e-3, Is=1e-6):
    w = 0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2)))
    if dIdt_richards(t, z) > 0:
        return vd_max*w
    else:
        return vr_max*w/(1+delta*P)

def v_richardsII(t, z, P, vd_max, vr_max, delta, Zmax, gamma=0.05, b=10):  # nimp
    w = 1 / (gamma * I_richards(t, z)) * dIdt_richards(t, z)
    if (w > vd_max):
        h = vd_max
    elif (w < -vr_max):
        h = vr_max
    else:
        h = w
    if 0 <= z <= Zmax - b:
        gd = 1
    else:
        gd = np.cos(np.pi / (2 * b) * (z - Zmax / b))
    if b <= z <= Zmax:
        gu = 1
    else:
        gu = np.sin(np.pi / (2 * b) * z)
    g = (0 <= t % 24 < 12) * gd / (1 + delta * P) + (12 <= t % 24 < 24) * gu
    W = h * g
    return W


def vc(t, z, P, c):
    th = t%24
    if th <= 12:
        return c
    else:
        return -c


def v_affine(t, z, P, a):
    th = t % 24
    return -a * th + a * 12


def v(t, z, P):
    th = t % 24
    if 4.5 <= th <= 6.5:
        return 72
    elif 16.5 <= th <= 18.5:
        return -72
    else:
        return 0


def v1(t, z, P, vd_max, vr_max, th, delta):
    if dIdt_richards(t, z) > th:
        return vd_max*(1-th/dIdt_richards(t, z))
    elif 0 <= dIdt_richards(t, z):
        return 0
    else:
        return vr_max * dIdt_richards(t, z) / (1 + delta*P)