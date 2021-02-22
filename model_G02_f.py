import numpy as np


def I0_richards(t, Is=1e-6, eps=1e-3):
    return Is + (1 - Is) / 2 * (
            1 + np.sin(np.pi / 12 * (t - 6)) + np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2) - np.sqrt(eps + 1))


def dIdt_richards(t, z, gamma=0.05, Is=1e-6, eps=1e-3):
    return np.exp(-gamma * z) * (0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2))))


def dI0dt_richards(t, Is=1e-6, eps=1e-3):
    return 0.5 * (1 - Is) * (
            np.pi / 12 * np.cos(np.pi / 12 * (t - 6)) * (1 + 1 / np.sqrt(eps + np.sin(np.pi / 12 * (t - 6)) ** 2)))


def I_richards(t, z, gamma=0.05):
    return I0_richards(t) * np.exp(- gamma * z)


def r(t, z, r_max, K_I, light):
    return r_max * light(t, z) / (K_I + light(t, z))


def v_richards_seuil(t, z, R, dt, tdI, v):
    """
    Migration speed
    :param t: time (h)
    :param z: depth (m)
    :param R: resource concentration
    :param tdI: treshold for the relative rate of light change under which the speed is null
    :param dt: time step
    :param v: constant speed (m/h)
    :return: migration speed at a given time and depth, can be positive (descent) or negative (ascent)
    """
    dI_r = 1/dt * (I_richards(t+dt, z)/I_richards(t, z)-1)
    if np.abs(dI_r) < tdI:
        return 0
    else:
        if dI_r > 0:
            return v
        else:
            return -v



def model_AN_RGCD(R0, C0, dz, dt, light, speed, args, K_I, r_max, alpha, beta, K, e, mu, d):
    """
    Numerical scheme from Andersen and Nival, 1991 with Cell Quota Model applied to a resource - migrating consumer
    1D model with the production of detritus
    :param R0: initialisation of the resource
    :param C0: initialisation of the consumer
    :param dz: depth step
    :param dt: time step
    :param light: light intensity as a function of time and depth
    :param speed: migration speed (positive downward, negative upward) with arguments (t, z, P, args)
    :param args: arguments for the speed function after (t, z, P), can be empty then put ()
    :param K_I: half saturation constant for the resource functional response to light
    :param r_max: maximum resource growth rate
    :param alpha: feeding rate
    :param beta: handling time
    :param K: carrying capacity for the resource
    :param e: assimilation coefficient/ gross growth efficiency
    :param mu: linear mortality rate for the consumer
    :param d: digestion rate
    :return: speed, resource, gut content, energy, consumer and detritus concentration in time and space
    """
    # initialisation
    tt, zz = np.shape(C0)
    C = np.copy(C0)
    R = np.copy(R0)
    G = np.zeros((tt, zz))
    D = np.zeros((tt, zz))
    w = np.zeros((tt, zz))
    for t in range(tt - 1):
        if (t * dt) % 24 == 0:
            print(t * dt)
        for i in range(zz):
            w[t, i] = speed(t * dt, i * dz, R[t, i], *args)
            # CFL condition
            if np.abs(w[t, i]) * dt / dz >= 1:
                print('CFL not satisfied : w*dt/dz=' + str(w[t, i]) + '*' + str(dt) + '/' + str(dz))
                return None

            R[t + 1, i] = R[t, i] + dt * (
                    r(t * dt, i * dz, r_max, K_I, light) * R[t, i] * (1 - R[t, i] / K)) - dt * alpha * R[t, i] * C[
                              t, i] / (1 + beta * R[t, i])
            D[t + 1, i] = D[t, i] + dt * d * (1 - e) * G[t, i] * C[t, i]  # + dt * mu * C[t, i]

        for i in range(1, zz - 1):
            G[t + 1, i] = G[t, i] + dt / dz * (
                    (w[t, i - 1] > 0) * w[t, i - 1] * G[t, i - 1] - np.abs(w[t, i]) * G[t, i] - 1 * (w[t, i + 1] < 0) *
                    w[t, i + 1] * G[t, i + 1]) + dt * alpha * R[t, i] * C[t, i] / (1 + beta * R[t, i]) - dt * d * G[t, i]
            C[t + 1, i] = C[t, i] + dt / dz * (
                        (w[t, i - 1] > 0) * w[t, i - 1] * C[t, i - 1] - np.abs(w[t, i]) * C[t, i] - 1 * (
                        w[t, i + 1] < 0) * w[t, i + 1] * C[t, i + 1]) + dt * d * e * G[t, i] - dt * mu * C[t, i]

        C[t + 1, 0] = C[t, 0] + dt / dz * (
                    -w[t, 0] * (w[t, 0] > 0) * C[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * C[t, 1]) + dt * e * d * G[t, 0] - dt * mu * C[t, 0]
        C[t + 1, -1] = C[t, -1] + dt / dz * (
                    (w[t, -2] > 0) * w[t, -2] * C[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * C[t, -1]) + dt * e * d * G[t, -1] - dt * mu * C[t, -1]

        G[t + 1, 0] = G[t, 0] + dt / dz * (
                -w[t, 0] * (w[t, 0] > 0) * G[t, 0] - 1 * (w[t, 1] < 0) * w[t, 1] * G[t, 1]) + dt * alpha * R[t, 0] / (
                              1 + beta * R[t, 0]) - dt * d * G[t, 0]
        G[t + 1, -1] = G[t, -1] + dt / dz * (
                (w[t, -2] > 0) * w[t, -2] * G[t, -2] - np.abs(w[t, -1]) * (w[t, -1] < 0) * G[t, -1]) + dt * alpha * \
                       R[t, -1] / (1 + beta * R[t, -1]) - dt * d * G[t, -1]

    return w, R, G, C, D