import numpy as np
import os

## Test the vectorized code:

def test():
    
    
    
    # row vector
    fs_7 = np.array([2, 1, 0, 0, 0, 0]).T
    
    # column vector
    b_7 = np.array([0, 0, 0, 0, ])


def free_surface_J_vect(**params):
    pass


def free_surface_J(**params):
    """
    Final equation in Method section. 
    """
    # vars to define, al, au, cp0, r, Du, Au, Bu, Cu
    dl = params["dl"]
    du = params["du"]

    rhou = params["rhou"]
    rhol = params["rhol"]

    mul = params["mul"]
    muu = params["muu"]

    g = params["g"]

    U0 = params["U0"]

    n = dl / du
    m = mul /muu
    r = rhol / rhou

    gamma = (n ** 3) * r / (m ** 2)
    Ku = params["Ku"]
    Kl = Ku / gamma

    Al = -(Kl * (du ** 2)) / (2 * U0 * mul)
    Au = m * gamma * Al

    al = (Ku * (du ** 2)) / (U0 * m * muu)
    au = m * al

    b = ((du ** 2) * n * Kl) * (gamma + .5 * n) / (U0 * mul)

    Fr2 = (rhol - rhou) * (g * du) / (rhol) / (U0 ** 2)
    Fr2 = np.sqrt(Fr2)

    Bu = - (2 * (n ** 3) + 3 * (n ** 2) + 2 * m) / ((n ** 2) * (2 * n + 3))
    Bl = 3 * (n + 2) / (n * (2 * n + 3))

    Cl = 3 / ((n ** 2) * (2 * n + 3))
    Cu = m * Cl

    Dl = -1 / ((n ** 2) * (2 * n + 3))
    Du = m * Dl

    cp0 = -((n ** 2)*(2 * n + 3) * (al - au)) / (2 * ((n ** 3) + (3 * n ** 2) + (3 * n) + m))

    hu = (Au * Du) / 210 + (au * Du) / 60 + (au * Cu - 3 * cp0 * Du - Au * Bu) / 60 + (cp0 * Cu + Au) / 12

    dhu2 = (Au * Du) / 5 + (au * Du) / 2 + (au * Cu - 3 * cp0 * Du - Au * Bu) / 3 + (cp0 * Cu + Au)

    hl = (Al * Dl) / 210 + (al * Dl) / 60 + (al * Cl - 3 * cp0 * Dl - Al * Bl) / 60 + (cp0 * Cl + Al) / 12

    hln = (Al * Dl) * (n ** 7) / 210 + (al * Dl) * (n ** 6) / 60 + (al * Cl - 3 * cp0 * Dl - Al * Bl) * (n ** 5) / 60 + (cp0 * Cl + Al) * (n ** 4)/ 12

    dhln = (Al * Dl)  * (n ** 6)/ 30 + (al * Dl) * (n ** 5)/ 10 + (al * Cl - 3 * cp0 * Dl - Al * Bl) * (n ** 4)/ 12 + (cp0 * Cl + Al) * (n ** 3)/ 3

    L = (n ** 3) * (Fr2 + cp0 * ((au + Bu*cp0) * (1 - r)))
    LL = 12 * cp0 * hl * r
    LLL = -cp0 * (n ** 3) * dhu2
    IL = 6 * cp0 * n * r * dhln
    G = 2 * cp0 * n * r * (3 * hln + dhln * n)
    GG = 3 * cp0 * (n ** 2) * (-dhu2 + 6 * hu + 4 * hu * n)

    J1 = cp0 / ((2 * (n ** 2) * m) * (2 * n + 3) * (al - au))
    J2 = L * ((4 * m - n) + LL * (m - n) + IL * (2 * m - n) + GG * m - G * n)

    J = J1 * J2
    return J
