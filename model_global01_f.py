import numpy as np


def model_AN_RC(P0, Z0, E0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu, rho, Em):
    tt, zz = np.shape(Z0)
    Z = np.copy(Z0)
    P = np.copy(P0)
    E = np.copy(E0)
    D = np.zeros((tt, zz))
    w = np.zeros((tt, zz))
    for t in range(tt - 1):
        if (t * dt)%24 == 0:
            print(t * dt)
        for i in range(zz):
            w[t, i] = speed(t * dt, i * dz, P[t, i], *args)
            # CFL condition
            if np.abs(w[t, i]) * dt / dz >= 1:
                print('CFL not satisfied : w*dt/dz=' + str(w[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            P[t + 1, i] = P[t, i] + dt * (
                        R(t * dt, i * dz, r_max, K_I, light) * P[t, i] * (1 - P[t, i] / K)) - dt * alpha * P[t, i] * Z[
                              t, i] / (1 + beta * P[t, i])
            D[t + 1, i] = D[t, i] + dt * (1 - e) * alpha * P[t, i] * Z[t, i] / (1 + beta * P[t, i]) + dt * mu * Z[t, i]

        for i in range(1, zz - 1):
            E[t + 1, i] = E[t, i] + dt / dz * (
                        (w[t, i - 1] > 0) * w[t, i - 1] * E[t, i - 1] - np.abs(w[t, i]) * E[t, i] - 1 * (
                        w[t, i + 1] < 0) * w[t, i + 1] * E[t, i + 1]) + dt * e * alpha * P[t, i] / (1 + beta * P[t, i]) - dt * rho * (
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
                    -w[t, 0] * (w[t, 0] > 0) * Z[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * Z[t, 1]) + dt * rho * \
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
                    -w[t, 0] * (w[t, 0] > 0) * E[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * E[t, 1]) + dt * e * alpha * P[t, 0] / (1 + beta * P[t, 0]) - dt * rho * (
                    E[t, 0] - Em)
        E[t + 1, -1] = E[t, -1] + dt / dz * (
                (w[t, -2] > 0) * w[t, -2] * E[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * E[t, -1]) + dt * e * alpha * P[t, -1] / (1 + beta * P[t, -1]) - dt * rho * (
                    E[t, -1] - Em)

    return w, P, E, Z, D


def I0_richards(t, Is=1e-6, eps=1e-4):
    return Is + (1 - Is) / 2 * (
            1 + np.sin(np.pi / 12 * (t - 6)) + np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2) - np.sqrt(eps + 1))


def dIdt_richards(t, z, gamma=0.05, Is=1e-6, eps=1e-4):
    return np.exp(-gamma * z) * (0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2))))


def dI0dt_richards(t, Is=1e-6, eps=1e-4):
    return 0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2)))


def I_richards(t, z, gamma=0.05):
    return I0_richards(t) * np.exp(- gamma * z)


def R(t, z, r_max, K_I, light):
    return r_max * light(t, z) / (K_I + light(t, z))


def v_richards_dI0dt(t, z, P, vd_max, vr_max, delta, tdI, td, tl, dl):
    """
    Migration speed
    :param t: time (h)
    :param z: depth (m)
    :param P: phytoplankton concentration
    :param vd_max: coefficient for descent speed
    :param vr_max: coefficient for ascent speed
    :param delta: handling time for feeding-dependent speed (during the ascent)
    :param tdI: treshold for the rate of surface light change under which the speed is null
    :param td: treshold in depth for lowering the ascent speed
    :param tl: treshold in absolute light for lowering the descent speed
    :param dl: such that the descent speed is null under tl - dl (absolute light)
    :return: numerical value for the migration time at a given time and depth, can be positive (descent) or negative (descent)
    """
    w = dI0dt_richards(t)
    if np.abs(w) < tdI:
        return 0
    elif w > 0:
        if tl - dl <= I_richards(t, z) <= tl:
            return np.cos(np.pi/2*(1-(I_richards(t, z)-(tl - dl))/dl))*vd_max*w
        elif I_richards(t, z) < tl - dl:
            return 0
        else:
            return vd_max*w
    else:
        if z == 0:
            return 0
        elif 0 < z < td:
            return np.exp(1-td/z)*vr_max*w/(1+delta*P)
        else:
            return vr_max*w/(1+delta*P)



