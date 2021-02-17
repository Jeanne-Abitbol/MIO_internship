import numpy as np


def model(func, args):
    return func(*args)


def norm(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def I0(t):
    th = t % 24
    if 4 <= th <= 20:
        return 750 * (1 + np.sin(2 * np.pi * th / 16 - np.pi))
    else:
        return t * 0


def I0_poggiale(t, alpha=1):
    th = t % 24
    if 4 < th < 20:
        return 1500 * np.exp((-alpha * (th - 12) ** 2) / ((th - 12) ** 2 - 64) ** 2)
    else:
        return 0


def I0_kamykowski(t):
    return 0.5 * (np.sin(np.pi / 12 * (t - 6)) + np.abs(np.sin(np.pi / 12 * (t - 6))))


def I(t, z, gamma=0.05):
    return I0(t) * np.exp(-gamma * z)


def I0_richards(t, Is=1e-6, eps=1e-4):
    return Is + (1 - Is) / 2 * (
            1 + np.sin(np.pi / 12 * (t - 6)) + np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2) - np.sqrt(eps + 1))


def dIdt_richards(t, z, gamma=0.05, Is=1e-6, eps=1e-4):
    return np.exp(-gamma * z) * (0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2))))


def dI0dt_richards(t, Is=1e-6, eps=1e-4):
    return (0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2))))


def I_richards(t, z, gamma=0.05):
    return I0_richards(t) * np.exp(- gamma * z)


def R(t, z, r_max, K_I, light):
    return r_max * light(t, z) / (K_I + light(t, z))


def model_AN_migrationonly(Z0, P0, dz, dt, speed, args):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    w = np.zeros((tt, zz))
    for t in range(tt):
        for i in range(zz):
            w[t, i] = speed(t * dt, i * dz, P[i], *args)
            # CFL condition
            if np.abs(w[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(w[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None
    for t in range(tt - 1):
        for i in range(1, zz - 1):
            Z[t + 1, i] = Z[t, i] + dt / dz * (
                    (w[t, i - 1] > 0) * w[t, i - 1] * Z[t, i - 1] - np.abs(w[t, i]) * Z[t, i] - 1 * (
                    w[t, i + 1] < 0) * w[t, i + 1] * Z[t, i + 1])
        Z[t + 1, 0] = Z[t, 0] + dt / dz * (
                -w[t, 0] * (w[t, 0] > 0) * Z[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * Z[t, 1])
        Z[t + 1, -1] = Z[t, -1] + dt / dz * (
                (w[t, -2] > 0) * w[t, -2] * Z[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * Z[t, -1])
    return w, Z


def model_sweby_migrationonly(Z0, P0, dz, dt, speed, args):
    """
    dB/dt = -wdB/dz
    :param Z0: 2D array
    :param P0: 1D array (distribution in the watercolumn constant in time)
    :param dz: depth step
    :param dt: time step
    :param speed:
    :param args:
    :return:
    """
    tt, zz = np.shape(Z0)
    # print(tt, zz)
    imax = zz - 1
    Z = np.copy(Z0)
    P = np.copy(P0)

    s = np.zeros((tt, zz))

    for t in range(tt - 1):
        nu = np.zeros(zz)
        r = np.zeros(zz)
        nu[0] = 1
        th = (t * dt) % 24
        if th == np.int(th):
            print(th)

        for i in range(zz):

            # computation of the migration speed
            s[t, i] = speed(t * dt, i * dz, P[i], *args)

            # CFL condition
            if np.abs(s[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

        for i in range(1, zz - 1):

            # the descent
            if s[t, i] >= 0:
                if Z[t, i + 1] - Z[t, i] != 0:
                    nu[i] = 1 - dt / dz * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i]) / (Z[t, i + 1] - Z[t, i])
                else:
                    nu[i] = 1 - dt / dz * (s[t, i + 1] - s[t, i])
                if (nu[i] * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i])) == 0:
                    r[i] = 1
                else:
                    r[i] = (nu[i - 1] * (s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1])) / (
                            nu[i] * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i]))  # day
                r[r < 0] = 0
                r[r > 1] = 1
                Z[t + 1, i] = Z[t, i] - dt / dz * (s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1]) + dt / (2 * dz) * (
                        r[i] * nu[i] * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i]) - r[i - 1] * nu[i - 1] * (
                        s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1]))

            # the ascension
            else:
                l = imax - i  # inversion of the axis as the scheme involves positive advective velocity
                if Z[t, l + 1] - Z[t, l] != 0:
                    nu[l] = 1 - dt / dz * ((np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l]) * Z[t, l]) / (
                            Z[t, l + 1] - Z[t, l]))
                else:
                    nu[l] = 1 - dt / dz * (np.abs(s[t, l + 1]) - np.abs(s[t, l]))
                if (nu[l] * (np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l]) * Z[t, l])) == 0:
                    r[l] = 1
                else:
                    r[l] = (nu[l - 1] * (np.abs(s[t, l]) * Z[t, l] - np.abs(s[t, l - 1]) * Z[t, l - 1])) / (
                            nu[l] * (np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l]) * Z[t, l]))
                r[r < 0] = 0
                r[r > 1] = 1
                Z[t + 1, l] = Z[t, l] - dt / dz * (
                        np.abs(s[t, l]) * Z[t, l] - np.abs(s[t, l + 1]) * Z[t, l + 1]) + dt / (2 * dz) * (
                                      r[l] * nu[l] * (
                                      np.abs(s[t, l - 1]) * Z[t, l - 1] - np.abs(s[t, l]) * Z[t, l]) - r[
                                          l + 1] * nu[l + 1] * (np.abs(s[t, l]) * Z[t, l] - np.abs(s[t, l + 1]) * Z[
                                  t, l + 1]))
        Z[t + 1, 0] = Z[t + 1, 2]
        Z[t + 1, imax] = Z[t + 1, imax - 2]

    return s, Z


def model_sweby(Z0, P0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu):
    """

    :param Z0: 2D array
    :param P0: 2D array
    :param dz: depth step
    :param dt: time step
    :param light: function of light with arguments t, z, and gamma as optional (fixed)
    :param speed: function for migration speed
    :param args: arguments following t, z, P, for the speed function
    :param K_I: half-saturation constant for phytoplankton functional response to light
    :param r_max: maximum phytoplankton growth rate
    :param alpha: maximum feeding rate
    :param beta: handling time for zooplankton feeding on phytoplankton
    :param K: capacitance (nutrients) for phytoplankton logistic growth
    :param e: assimilation coefficient
    :param mu: linear mortality rate for zooplankton
    :return: s, P, Z
    """
    tt, zz = np.shape(Z0)
    # print(tt, zz)
    imax = zz - 1
    Z = np.copy(Z0)
    P = np.copy(P0)

    s = np.zeros((tt, zz))

    for t in range(tt - 1):
        nu = np.zeros(zz)
        r = np.zeros(zz)
        nu[0] = 1
        th = (t * dt) % 24
        if th == np.int(th):
            print(th)

        for i in range(zz):
            # phytoplankton
            P[t + 1, i] = P[t, i] + dt * (
                    R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])

            # computation of the migration speed
            s[t, i] = speed(t * dt, i * dz, P[t, i], *args)

            # CFL condition
            if np.abs(s[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

        for i in range(1, zz - 1):

            # the descent
            if s[t, i] >= 0:
                if Z[t, i + 1] - Z[t, i] != 0:
                    nu[i] = 1 - dt / dz * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i]) / (Z[t, i + 1] - Z[t, i])
                else:
                    nu[i] = 1 - dt / dz * (s[t, i + 1] - s[t, i])
                if (nu[i] * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i])) == 0:
                    r[i] = 1
                else:
                    r[i] = (nu[i - 1] * (s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1])) / (
                            nu[i] * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i]))  # day
                r[r < 0] = 0
                r[r > 1] = 1
                Z[t + 1, i] = Z[t, i] - dt / dz * (s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1]) + dt / (2 * dz) * (
                        r[i] * nu[i] * (s[t, i + 1] * Z[t, i + 1] - s[t, i] * Z[t, i]) - r[i - 1] * nu[i - 1] * (
                        s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1])) + dt * e * alpha * P[t, i] * Z[t, i] / (
                                      1 + beta * P[t, i]) - dt * mu * Z[t, i]

            # the ascension
            else:
                l = imax - i  # inversion of the axis as the scheme involves positive advective velocity
                if Z[t, l + 1] - Z[t, l] != 0:
                    nu[l] = 1 - dt / dz * ((np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l]) * Z[t, l]) / (
                            Z[t, l + 1] - Z[t, l]))
                else:
                    nu[l] = 1 - dt / dz * (np.abs(s[t, l + 1]) - np.abs(s[t, l]))
                if (nu[l] * (np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l]) * Z[t, l])) == 0:
                    r[l] = 1
                else:
                    r[l] = (nu[l - 1] * (np.abs(s[t, l]) * Z[t, l] - np.abs(s[t, l - 1]) * Z[t, l - 1])) / (
                            nu[l] * (np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l]) * Z[t, l]))
                r[r < 0] = 0
                r[r > 1] = 1
                Z[t + 1, l] = Z[t, l] - dt / dz * (
                        np.abs(s[t, l]) * Z[t, l] - np.abs(s[t, l + 1]) * Z[t, l + 1]) + dt / (2 * dz) * (
                                      r[l] * nu[l] * (
                                      np.abs(s[t, l - 1]) * Z[t, l - 1] - np.abs(s[t, l]) * Z[t, l]) - r[
                                          l + 1] * nu[l + 1] * (np.abs(s[t, l]) * Z[t, l] - np.abs(s[t, l + 1]) * Z[
                                  t, l + 1])) + dt * e * alpha * P[t, l] * Z[t, l] / (1 + beta * P[t, l]) - dt * mu * Z[
                                  t, l]

        Z[t + 1, 0] = Z[t + 1, 2]
        Z[t + 1, imax] = Z[t + 1, imax - 2]

    return s, P, Z


def model_AN(Z0, P0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    w = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        for i in range(zz):
            w[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            # CFL condition
            if np.abs(w[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(w[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            P[t + 1, i] = P[t, i] + dt * (
                        R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])

        for i in range(1, zz - 1):
            Z[t + 1, i] = Z[t, i] + dt / dz * (
                    (w[t, i - 1] > 0) * w[t, i - 1] * Z[t, i - 1] - np.abs(w[t, i]) * Z[t, i] - 1 * (
                    w[t, i + 1] < 0) * w[t, i + 1] * Z[t, i + 1]) + dt * e * alpha * P[t, i] * Z[t, i] / (
                                  1 + beta * P[t, i]) - dt * mu * Z[t, i]
        Z[t + 1, 0] = Z[t, 0] + dt / dz * (-np.abs(w[t, 0]) * (w[t, 0] > 0) * Z[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * Z[
            t, 1]) + dt * e * alpha * P[t, 0] * Z[t, 0] / (1 + beta * P[t, 0]) - dt * mu * Z[t, 0]
        Z[t + 1, -1] = Z[t, -1] + dt / dz * (
                (w[t, -2] > 0) * w[t, -2] * Z[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * Z[
            t, -1]) + dt * e * alpha * P[t, -1] * Z[t, -1] / (1 + beta * P[t, -1]) - dt * mu * Z[t, -1]

    return w, P, Z


def model_AN_RC(P0, Z0, E0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu, rho, Em):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    E = np.copy(E0)
    D = np.zeros((tt, zz))
    w = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        for i in range(zz):
            w[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            # CFL condition
            if np.abs(w[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(w[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            P[t + 1, i] = P[t, i] + dt * (
                        R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])
            D[t + 1, i] = D[t, i] + dt * (1 - e) * alpha * P[t, i] * Z[t, i] / (1 + beta * P[t, i]) + dt * mu * Z[t, i]

        for i in range(1, zz - 1):
            E[t + 1, i] = E[t, i] + dt / dz * (
                        (w[t, i - 1] > 0) * w[t, i - 1] * E[t, i - 1] - np.abs(w[t, i]) * E[t, i] - 1 * (
                        w[t, i + 1] < 0) * w[t, i + 1] * E[t, i + 1]) + dt * e * alpha * P[t, i] * Z[t, i] / (1 + beta * P[t, i]) - dt * rho * (
                    E[t, i] - Em)
            if E[t, i] > 0:
                Z[t + 1, i] = Z[t, i] + dt / dz * (
                        (w[t, i - 1] > 0) * w[t, i - 1] * Z[t, i - 1] - np.abs(w[t, i]) * Z[t, i] - 1 * (
                        w[t, i + 1] < 0) * w[t, i + 1] * Z[t, i + 1]) + dt * rho * Z[t, i] * (
                                      1 - Em / E[t, i]) - dt * mu * Z[t, i]
            else:
                Z[t + 1, i] = 0
        if E[t, 0] > 0:
            Z[t + 1, 0] = Z[t, 0] + dt / dz * (
                    -np.abs(w[t, 0]) * (w[t, 0] > 0) * Z[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * Z[t, 1]) + dt * rho * \
                          Z[t, 0] * (1 - Em / E[t, 0]) - dt * mu * Z[t, 0]
        else:
            Z[t + 1, 0] = 0
        if E[t, -1] > 0:
            Z[t + 1, -1] = Z[t, -1] + dt / dz * (
                    (w[t, -2] > 0) * w[t, -2] * Z[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * Z[t, -1]) + dt * rho * Z[
                               t, 0] * (1 - Em / E[t, -1]) - dt * mu * Z[t, -1]
        else:
            Z[t + 1, -1] = 0
        E[t + 1, 0] = E[t, 0] + dt / dz * (
                -np.abs(w[t, 0]) * (w[t, 0] > 0) * E[t, 0] - 1 * (w[t, 1] < 0) * E[t, 1] * Z[t, 1]) + dt * e * alpha * P[t, 0] * Z[t, 0] / (1 + beta * P[t, 0]) - dt * rho * (
                    E[t, 0] - Em)
        E[t + 1, -1] = E[t, -1] + dt / dz * (
                (w[t, -2] > 0) * w[t, -2] * E[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * E[t, -1]) + dt * e * alpha * P[t, -1] * Z[t, -1] / (1 + beta * P[t, -1]) - dt * rho * (
                    E[t, -1] - Em)

    return w, P, Z, D


def model_AN_RC_withoutE(P0, Z0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    D = np.zeros((tt, zz))
    w = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        for i in range(zz):
            w[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            # CFL condition
            if np.abs(w[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(w[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            P[t + 1, i] = P[t, i] + dt * (
                        R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])
            D[t + 1, i] = D[t, i] + dt * (1 - e) * alpha * P[t, i] * Z[t, i] / (1 + beta * P[t, i]) + dt * mu * Z[t, i]

        for i in range(1, zz - 1):
            Z[t + 1, i] = Z[t, i] + dt / dz * (
                    (w[t, i - 1] > 0) * w[t, i - 1] * Z[t, i - 1] - np.abs(w[t, i]) * Z[t, i] - 1 * (
                    w[t, i + 1] < 0) * w[t, i + 1] * Z[t, i + 1]) + dt * e * alpha * P[t, i] * Z[t, i] / (
                                      1 + beta * P[t, i]) - dt * mu * Z[t, i]
        Z[t + 1, 0] = Z[t, 0] + dt / dz * (
                -np.abs(w[t, 0]) * (w[t, 0] > 0) * Z[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * Z[t, 1]) + dt * e * alpha * \
                      P[t, 0] * Z[t, 0] / (1 + beta * P[t, 0]) - dt * mu * Z[t, 0]
        Z[t + 1, -1] = Z[t, -1] + dt / dz * (
                (w[t, -2] > 0) * w[t, -2] * Z[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * Z[t, -1]) + dt * e * alpha * \
                       P[t, -1] * Z[t, -1] / (1 + beta * P[t, -1]) - dt * mu * Z[t, -1]

    return w, P, Z, D


def leapfrog(Z0, P0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    s = np.zeros((tt, zz))
    t = 0
    for i in range(zz):
        s[t, i] = speed(t * dt, i * dz, P[t, i], *args)
        # CFL condition
        if np.abs(s[t, i]) * dt / dz > 1:
            print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
            return None

        P[t + 1, i] = P[t, i] + dt * (
                R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                          t, i] / (1 + beta * P[t, i])
    for i in range(1, zz - 1):
        Z[t + 1, i] = Z[t, i] - dt / dz * (s[t, i] * Z[t, i] - s[t, i - 1] * Z[t, i - 1]) + dt * e * alpha * \
                      Z[t, i] * P[t, i] / (1 + beta * P[t, i]) - dt * mu * Z[t, i]
    Z[t + 1, 0] = Z[t + 1, 2]
    Z[t + 1, zz - 1] = Z[t + 1, zz - 3]

    for t in range(1, tt - 1):
        if t * dt == np.int(t * dt):
            print(t * dt)
        for i in range(zz):
            s[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            # CFL condition
            if np.abs(s[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            P[t + 1, i] = P[t, i] + dt * (
                    R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])
        for i in range(1, zz - 1):
            if s[t, i] >= 0:
                Z[t + 1, i] = Z[t - 1, i] - dt / dz * (
                        s[t, i + 1] * Z[t, i + 1] - s[t, i - 1] * Z[t, i - 1]) + dt * e * alpha * \
                              Z[t, i] * P[t, i] / (1 + beta * P[t, i]) - dt * mu * Z[t, i]
            else:
                l = zz - 1 - i
                Z[t + 1, l] = Z[t - 1, l] - dt / dz * (
                        np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l - 1]) * Z[t, l - 1]) + dt * e * alpha * \
                              Z[t, l] * P[t, l] / (1 + beta * P[t, l]) - dt * mu * Z[t, l]
        Z[t + 1, 0] = Z[t + 1, 2]
        Z[t + 1, zz - 1] = Z[t + 1, zz - 3]

    return s, P, Z


def leapfrog_migrationonly(Z0, P0, dz, dt, speed, args):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    s = np.zeros((tt, zz))
    t = 0
    for i in range(zz):
        s[t, i] = speed(t * dt, i * dz, P[i], *args)
        # CFL condition
        if np.abs(s[t, i]) * dt / dz > 1:
            print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
            return None

    Z[t + 1, 0] = Z[t, 0] + dt / dz * (
            -s[t, 0] * (s[t, 0] > 0) * Z[t, 0] - 1 * (s[t, 1] < 0) * s[t, 1] * Z[t, 1])
    Z[t + 1, -1] = Z[t, -1] + dt / dz * (
            (s[t, -2] > 0) * s[t, -2] * Z[t, -2] - np.abs(s[t, -1]) * (s[t, -1] < 0) * Z[t, -1])

    for i in range(1, zz - 1):
        Z[t + 1, i] = Z[t, i] + dt / dz * (
                (s[t, i - 1] > 0) * s[t, i - 1] * Z[t, i - 1] - np.abs(s[t, i]) * Z[t, i] - 1 * (
                s[t, i + 1] < 0) * s[t, i + 1] * Z[t, i + 1])

    for t in range(1, tt - 1):
        if (t * dt) % 24 == 0:
            print((t * dt)/24)
        for i in range(zz):
            s[t, i] = speed(t * dt, i * dz, P[i], *args)
            # CFL condition
            if np.abs(s[t, i]) * dt / dz > 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

        for i in range(1, zz - 1):
#            if s[t, i] >= 0:
            Z[t + 1, i] = Z[t - 1, i] - dt / dz * (s[t, i + 1] * Z[t, i + 1] - s[t, i - 1] * Z[t, i - 1])
#            else:
#                l = zz - 1 - i
#                Z[t + 1, l] = Z[t - 1, l] - dt / dz * (np.abs(s[t, l + 1]) * Z[t, l + 1] - np.abs(s[t, l - 1]) * Z[t, l - 1])
        Z[t + 1, 0] = Z[t, 0] + dt / dz * (
                -s[t, 0] * (s[t, 0] > 0) * Z[t, 0] - 1 * (s[t, 1] < 0) * s[t, 1] * Z[t, 1])
        Z[t + 1, -1] = Z[t, -1] + dt / dz * (
                (s[t, -2] > 0) * s[t, -2] * Z[t, -2] - np.abs(s[t, -1]) * (s[t, -1] < 0) * Z[t, -1])

    return s, Z


def transport(Z0, P0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu, schema):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    s = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        a = np.zeros((zz))
        b = np.zeros((zz))
        c = np.zeros((zz))
        for i in range(zz):
            s[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            nu = s[t, i] * dt / dz
            # CFL condition
            if np.abs(nu) > 1:
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

            # phytoplankton
            P[t + 1, i] = P[t, i] + dt * (
                    R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])

        for i in range(1, zz - 1):
            Z[t + 1, i] = a[i] * Z[t, i + 1] + b[i] * Z[t, i] + c[i] * Z[t, i - 1] + dt * e * alpha * P[t, i] * Z[
                t, i] / (1 + beta * P[t, i]) - dt * mu * Z[t, i]
        Z[t + 1, 0] = Z[t + 1, 1]
        Z[t + 1, -1] = Z[t + 1, -2]

    return s, P, Z


def transport2(Z0, P0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    s = np.zeros((tt, zz))
    for t in range(tt - 1):
        # if t * dt == np.int(t * dt):
        # print(t * dt)
        for i in range(zz):
            s[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            nu = s[t, i] * dt / dz
            # CFL condition
            if np.abs(nu) > 1:
                print('CFL not satisfied : w*dt/dz=' + str(s[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            # phytoplankton
            P[t + 1, i] = P[t, i] + dt * (
                    R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])

        for i in range(1, zz - 1):
            Z[t + 1, i] = Z[t, i] + dt / dz * ((s[t, i] > 0) * s[t, i] * (Z[t, i] - Z[t, i - 1]) - (s[t, i] < 0) *
                                               s[t, i] * (Z[t, i + 1] - Z[t, i])) + dt * e * alpha * \
                          P[t, i] * Z[t, i] / (1 + beta * P[t, i]) - dt * mu * Z[t, i]
        Z[t + 1, 0] = Z[t + 1, 1]
        Z[t + 1, -1] = Z[t + 1, -2]

    return s, P, Z
