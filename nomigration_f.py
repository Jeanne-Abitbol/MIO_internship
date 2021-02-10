from model_functions import *


def PZ_nomigration(Y, t, z, light, K_I, r_max, alpha, beta, K, e, mu):
    P = Y[0]
    Z = Y[1]
    dP = R(t, z, r_max, K_I, light) * P * (1 - P / K) - alpha * P * Z / (1 + beta * P)
    dZ = Z * (e * alpha * P / (1 + beta * P) - mu)
    return dP, dZ

def PZ_nomigration_withD(Y, t, z, light, K_I, r_max, alpha, beta, K, e, mu):
    P = Y[0]
    Z = Y[1]
    dP = R(t, z, r_max, K_I, light) * P * (1 - P / K) - alpha * P * Z / (1 + beta * P)
    dZ = Z * (e * alpha * P / (1 + beta * P) - mu)
    dD = (1-e) * alpha * P * Z / (1 + beta * P) + mu*Z
    return dP, dZ, dD


# variables P (phytoplankton), Z (zooplankton), N (micronekton), D(detritus)
def PZNDdroop(Y, t, z, K, e, r_max, K_I, light, alpha_PZ, alpha_PN, alpha_ZN, beta_PZ, beta_PN, beta_ZN, m_P, m_Z, m_N,
              mu_N, rho_Z, rho_N, EZm, ENm):
    """

    :param Y: Vector of state variables P, EZ, Z, EN, N, D
    :param t: time
    :param z: depth
    :param K: environement capacitance for phytoplankton growth
    :param e: assimilation coefficient
    :param r_max: maximum phytoplankton growth rate
    :param K_I: half-saturation constant for phytoplankton light-dependent functional response
    :param light: function of the time and depth for light intensity
    :param alpha_PZ: maximum feeding rate of zooplankton over phytoplankton
    :param alpha_PN: maximum feeding rate of micronekton over phytoplankton
    :param alpha_ZN: maximum feeding rate for micronekton over zooplankton
    :param beta_PZ: handling time for zooplankton feeding over phytoplankton
    :param beta_PN: handling time for micronekton feeding over phytoplankton
    :param beta_ZN: handling time for micronekton over zooplankton
    :param m_P: maintenance rate for phytoplankton
    :param m_Z: maintenance rate for zooplankton
    :param m_N: maintenance rate for micronekton
    :param mu_N: mortality rate for micronekton
    :param rho_Z: coefficient for zooplankton growth rate
    :param rho_N: coefficient for micronekton growth rate
    :param EZm: basal metabolic rate
    :param ENm: basal metabolic rate
    :return:
    """
    P, EZ, Z, EN, N, D = Y
    dP = R(t, z, r_max, K_I, light) * P * (1 - P / K) - alpha_PZ * P * Z / (1 + beta_PZ * P) - alpha_PN * P * N / (
                1 + beta_PN * P) - m_P * P
    print(dP)
    dEZ = e * alpha_PZ * P / (1 + beta_PZ * P) - rho_Z * (EZ - EZm)
    dZ = rho_Z * Z * (1 - EZm / EZ) - alpha_ZN * Z * N / (1 + beta_ZN * Z) - m_Z * Z
    dEN = e * (alpha_PN * P / (1 + beta_PN * P) + alpha_ZN * Z / (1 + beta_ZN * Z)) - rho_N * (EN - ENm)
    dN = rho_N * N * (1 - ENm / EN) - (m_N + mu_N) * N
    dD = (1 - e) * (alpha_PZ * P * Z / (1 + beta_PZ * P) + alpha_PN * P * N / (1 + beta_PN * P) + alpha_ZN * Z * N / (
                1 + beta_ZN * Z)) + m_P * P + m_Z * Z + m_N * N
    return dP, dEZ, dZ, dEN, dN, dD

def RC_droop(Y, t, z, light, K_I, r_max, alpha, beta, K, e, mu, rho, Em):
    RP = Y[0]
    E = Y[1]
    C = Y[2]
    dRP = R(t, z, r_max, K_I, light)*RP*(1-RP/K)-alpha*RP*C/(1+beta*RP)
    dE = e*alpha*RP*C/(1+beta*RP) - rho*(E-Em)
    dC = rho*C*(1-Em/E) - mu*C
    return dRP, dE, dC


def RC_droop_withD(Y, t, z, light, K_I, r_max, alpha, beta, K, e, mu, rho, Em):
    RP = Y[0] # ressource primaire
    E = Y[1] # réserve d'énergie
    C = Y[2] # consommateur # detritus
    dRP = R(t, z, r_max, K_I, light)*RP*(1-RP/K)-alpha*RP*C/(1+beta*RP)
    dE = e*alpha*RP*C/(1+beta*RP) - rho*(E-Em)
    dC = rho*C*(1-Em/E) - mu*C
    dD = (1-e)*alpha*RP*C/(1+beta*RP)+mu*C
    return dRP, dE, dC, dD

def phyto(Y, t, z, light, K_I, r_max, K):
    return R(t, z, r_max, K_I, light) * Y * (1 - Y / K)

