import numpy as np
import os
import matplotlib.pyplot as plt
import math
from itertools import cycle

"""
Taken from Matlab function solutions that can be found within the 'Free Slip LSA' directory.
"""

def no_slip_J(**params):
    
    Ku_p, n, mu_u, m, U_p, n, du, F, r, Kl_p, m, mu_l = \
        (params["Ku_p"], params["n"], params["mu_u"], params["m"], params["U_p"], 
        params["n"], params["du"], params["F"], params["r"], params["Kl_p"], 
        params["m"], params["mu_l"])
    
    return -(560*F**2*U_p**2*mu_l**2*mu_u**2*n**14 + 56*Kl_p**2*du**4*m**4*mu_u**2*n**8 + 28*Kl_p**2*du**4*m**5*mu_u**2*n**7 + 7*Kl_p**2*du**4*m**6*mu_u**2*n**6 + 56*Kl_p**2*du**4*m**3*mu_u**2*n**10 + 252*Kl_p**2*du**4*m**4*mu_u**2*n**9 + 84*Kl_p**2*du**4*m**5*mu_u**2*n**8 - 7*Kl_p**2*du**4*m**7*mu_u**2*n**6 - 672*Kl_p**2*du**4*m**2*mu_u**2*n**12 - 1372*Kl_p**2*du**4*m**3*mu_u**2*n**11 - 854*Kl_p**2*du**4*m**4*mu_u**2*n**10 - 336*Kl_p**2*du**4*m**5*mu_u**2*n**9 - 140*Kl_p**2*du**4*m**6*mu_u**2*n**8 - 28*Kl_p**2*du**4*m**7*mu_u**2*n**7 - 1372*Kl_p**2*du**4*m**2*mu_u**2*n**13 - 1988*Kl_p**2*du**4*m**3*mu_u**2*n**12 - 224*Kl_p**2*du**4*m**4*mu_u**2*n**11 + 686*Kl_p**2*du**4*m**5*mu_u**2*n**10 + 84*Kl_p**2*du**4*m**6*mu_u**2*n**9 - 945*Kl_p**2*du**4*m**2*mu_u**2*n**14 - 336*Kl_p**2*du**4*m**3*mu_u**2*n**13 + 1708*Kl_p**2*du**4*m**4*mu_u**2*n**12 + 1596*Kl_p**2*du**4*m**5*mu_u**2*n**11 + 112*Kl_p**2*du**4*m**6*mu_u**2*n**10 - 224*Kl_p**2*du**4*m**2*mu_u**2*n**15 + 553*Kl_p**2*du**4*m**3*mu_u**2*n**14 + 1484*Kl_p**2*du**4*m**4*mu_u**2*n**13 + 952*Kl_p**2*du**4*m**5*mu_u**2*n**12 + 224*Kl_p**2*du**4*m**3*mu_u**2*n**15 + 392*Kl_p**2*du**4*m**4*mu_u**2*n**14 + 224*Kl_p**2*du**4*m**5*mu_u**2*n**13 + 32*Ku_p**2*du**4*m**4*mu_l**2*n**4 - 2*Ku_p**2*du**4*m**5*mu_l**2*n**3 + Ku_p**2*du**4*m**6*mu_l**2*n**2 + 72*Ku_p**2*du**4*m**3*mu_l**2*n**6 + 148*Ku_p**2*du**4*m**4*mu_l**2*n**5 - 34*Ku_p**2*du**4*m**5*mu_l**2*n**4 + 8*Ku_p**2*du**4*m**6*mu_l**2*n**3 - Ku_p**2*du**4*m**7*mu_l**2*n**2 - 216*Ku_p**2*du**4*m**2*mu_l**2*n**8 - 192*Ku_p**2*du**4*m**3*mu_l**2*n**7 - 154*Ku_p**2*du**4*m**4*mu_l**2*n**6 - 240*Ku_p**2*du**4*m**5*mu_l**2*n**5 + 2*Ku_p**2*du**4*m**6*mu_l**2*n**4 - 6*Ku_p**2*du**4*m**7*mu_l**2*n**3 - 292*Ku_p**2*du**4*m**2*mu_l**2*n**9 + 64*Ku_p**2*du**4*m**3*mu_l**2*n**8 + 48*Ku_p**2*du**4*m**4*mu_l**2*n**7 - 6*Ku_p**2*du**4*m**5*mu_l**2*n**6 + 92*Ku_p**2*du**4*m**6*mu_l**2*n**5 + 361*Ku_p**2*du**4*m**2*mu_l**2*n**10 + 1040*Ku_p**2*du**4*m**3*mu_l**2*n**9 + 464*Ku_p**2*du**4*m**4*mu_l**2*n**8 + 144*Ku_p**2*du**4*m**5*mu_l**2*n**7 + 88*Ku_p**2*du**4*m**6*mu_l**2*n**6 + 264*Ku_p**2*du**4*m**2*mu_l**2*n**11 + 199*Ku_p**2*du**4*m**3*mu_l**2*n**10 - 524*Ku_p**2*du**4*m**4*mu_l**2*n**9 - 312*Ku_p**2*du**4*m**5*mu_l**2*n**8 - 98*Ku_p**2*du**4*m**2*mu_l**2*n**12 - 378*Ku_p**2*du**4*m**3*mu_l**2*n**11 - 560*Ku_p**2*du**4*m**4*mu_l**2*n**10 - 224*Ku_p**2*du**4*m**5*mu_l**2*n**9 + 114*Ku_p**2*du**4*m*mu_l**2*n**11 + 98*Ku_p**2*du**4*m*mu_l**2*n**12 - 6*Kl_p**2*du**4*mu_u**2*n**17*r - Kl_p**2*du**4*mu_u**2*n**18*r - 28*Ku_p**2*du**4*mu_l**2*n**13*r - 7*Ku_p**2*du**4*mu_l**2*n**14*r + 88*Kl_p**2*du**4*m*mu_u**2*n**14*r + 92*Kl_p**2*du**4*m*mu_u**2*n**15*r + 2*Kl_p**2*du**4*m*mu_u**2*n**16*r + 8*Kl_p**2*du**4*m*mu_u**2*n**17*r + Kl_p**2*du**4*m*mu_u**2*n**18*r + 112*Ku_p**2*du**4*m*mu_l**2*n**10*r + 84*Ku_p**2*du**4*m*mu_l**2*n**11*r - 140*Ku_p**2*du**4*m*mu_l**2*n**12*r + 7*Ku_p**2*du**4*m*mu_l**2*n**14*r + 4480*F**2*U_p**2*m*mu_l**2*mu_u**2*n**11 + 6720*F**2*U_p**2*m*mu_l**2*mu_u**2*n**12 + 6160*F**2*U_p**2*m*mu_l**2*mu_u**2*n**13 - 224*Kl_p**2*du**4*m**2*mu_u**2*n**11*r - 560*Kl_p**2*du**4*m**3*mu_u**2*n**10*r - 378*Kl_p**2*du**4*m**4*mu_u**2*n**9*r - 98*Kl_p**2*du**4*m**5*mu_u**2*n**8*r - 312*Kl_p**2*du**4*m**2*mu_u**2*n**12*r - 524*Kl_p**2*du**4*m**3*mu_u**2*n**11*r + 199*Kl_p**2*du**4*m**4*mu_u**2*n**10*r + 264*Kl_p**2*du**4*m**5*mu_u**2*n**9*r + 98*Kl_p**2*du**4*m**6*mu_u**2*n**8*r + 144*Kl_p**2*du**4*m**2*mu_u**2*n**13*r + 464*Kl_p**2*du**4*m**3*mu_u**2*n**12*r + 1040*Kl_p**2*du**4*m**4*mu_u**2*n**11*r + 361*Kl_p**2*du**4*m**5*mu_u**2*n**10*r + 114*Kl_p**2*du**4*m**6*mu_u**2*n**9*r - 6*Kl_p**2*du**4*m**2*mu_u**2*n**14*r + 48*Kl_p**2*du**4*m**3*mu_u**2*n**13*r + 64*Kl_p**2*du**4*m**4*mu_u**2*n**12*r - 292*Kl_p**2*du**4*m**5*mu_u**2*n**11*r - 240*Kl_p**2*du**4*m**2*mu_u**2*n**15*r - 154*Kl_p**2*du**4*m**3*mu_u**2*n**14*r - 192*Kl_p**2*du**4*m**4*mu_u**2*n**13*r - 216*Kl_p**2*du**4*m**5*mu_u**2*n**12*r - 34*Kl_p**2*du**4*m**2*mu_u**2*n**16*r + 148*Kl_p**2*du**4*m**3*mu_u**2*n**15*r + 72*Kl_p**2*du**4*m**4*mu_u**2*n**14*r - 2*Kl_p**2*du**4*m**2*mu_u**2*n**17*r + 32*Kl_p**2*du**4*m**3*mu_u**2*n**16*r + 224*Ku_p**2*du**4*m**2*mu_l**2*n**7*r + 392*Ku_p**2*du**4*m**3*mu_l**2*n**6*r + 224*Ku_p**2*du**4*m**4*mu_l**2*n**5*r + 952*Ku_p**2*du**4*m**2*mu_l**2*n**8*r + 1484*Ku_p**2*du**4*m**3*mu_l**2*n**7*r + 553*Ku_p**2*du**4*m**4*mu_l**2*n**6*r - 224*Ku_p**2*du**4*m**5*mu_l**2*n**5*r + 1596*Ku_p**2*du**4*m**2*mu_l**2*n**9*r + 1708*Ku_p**2*du**4*m**3*mu_l**2*n**8*r - 336*Ku_p**2*du**4*m**4*mu_l**2*n**7*r - 945*Ku_p**2*du**4*m**5*mu_l**2*n**6*r + 686*Ku_p**2*du**4*m**2*mu_l**2*n**10*r - 224*Ku_p**2*du**4*m**3*mu_l**2*n**9*r - 1988*Ku_p**2*du**4*m**4*mu_l**2*n**8*r - 1372*Ku_p**2*du**4*m**5*mu_l**2*n**7*r - 336*Ku_p**2*du**4*m**2*mu_l**2*n**11*r - 854*Ku_p**2*du**4*m**3*mu_l**2*n**10*r - 1372*Ku_p**2*du**4*m**4*mu_l**2*n**9*r - 672*Ku_p**2*du**4*m**5*mu_l**2*n**8*r + 84*Ku_p**2*du**4*m**2*mu_l**2*n**12*r + 252*Ku_p**2*du**4*m**3*mu_l**2*n**11*r + 56*Ku_p**2*du**4*m**4*mu_l**2*n**10*r + 28*Ku_p**2*du**4*m**2*mu_l**2*n**13*r + 56*Ku_p**2*du**4*m**3*mu_l**2*n**12*r + 8960*F**2*U_p**2*m**2*mu_l**2*mu_u**2*n**8 + 31360*F**2*U_p**2*m**3*mu_l**2*mu_u**2*n**7 + 40880*F**2*U_p**2*m**4*mu_l**2*mu_u**2*n**6 + 24080*F**2*U_p**2*m**5*mu_l**2*mu_u**2*n**5 + 6160*F**2*U_p**2*m**6*mu_l**2*mu_u**2*n**4 + 560*F**2*U_p**2*m**7*mu_l**2*mu_u**2*n**3 + 26880*F**2*U_p**2*m**2*mu_l**2*mu_u**2*n**9 + 87360*F**2*U_p**2*m**3*mu_l**2*mu_u**2*n**8 + 100800*F**2*U_p**2*m**4*mu_l**2*mu_u**2*n**7 + 47040*F**2*U_p**2*m**5*mu_l**2*mu_u**2*n**6 + 6720*F**2*U_p**2*m**6*mu_l**2*mu_u**2*n**5 + 52640*F**2*U_p**2*m**2*mu_l**2*mu_u**2*n**10 + 135520*F**2*U_p**2*m**3*mu_l**2*mu_u**2*n**9 + 135520*F**2*U_p**2*m**4*mu_l**2*mu_u**2*n**8 + 52640*F**2*U_p**2*m**5*mu_l**2*mu_u**2*n**7 + 4480*F**2*U_p**2*m**6*mu_l**2*mu_u**2*n**6 + 47040*F**2*U_p**2*m**2*mu_l**2*mu_u**2*n**11 + 100800*F**2*U_p**2*m**3*mu_l**2*mu_u**2*n**10 + 87360*F**2*U_p**2*m**4*mu_l**2*mu_u**2*n**9 + 26880*F**2*U_p**2*m**5*mu_l**2*mu_u**2*n**8 + 24080*F**2*U_p**2*m**2*mu_l**2*mu_u**2*n**12 + 40880*F**2*U_p**2*m**3*mu_l**2*mu_u**2*n**11 + 31360*F**2*U_p**2*m**4*mu_l**2*mu_u**2*n**10 + 8960*F**2*U_p**2*m**5*mu_l**2*mu_u**2*n**9 - 88*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**6 - 26*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**5 - 8*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**4 - 128*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**8 - 400*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**7 - 50*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**6 - 8*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**5 + 8*Kl_p*Ku_p*du**4*m**7*mu_l*mu_u*n**4 + 888*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**10 + 1564*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**9 + 1008*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**8 + 576*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**7 + 138*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**6 + 34*Kl_p*Ku_p*du**4*m**7*mu_l*mu_u*n**5 + 1664*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**11 + 1924*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**10 + 176*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**9 - 680*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**8 - 176*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**7 + 584*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**12 - 704*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**11 - 2172*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**10 - 1740*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**9 - 200*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**8 - 40*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**13 - 752*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**12 - 960*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**11 - 640*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**10 + 98*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**14 + 154*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**13 + 168*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**12 - 114*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**13 - 98*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**14 + 34*Kl_p*Ku_p*du**4*mu_l*mu_u*n**15*r + 8*Kl_p*Ku_p*du**4*mu_l*mu_u*n**16*r - 200*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**12*r - 176*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**13*r + 138*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**14*r - 8*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**15*r - 8*Kl_p*Ku_p*du**4*m*mu_l*mu_u*n**16*r + 168*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**8*r + 154*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**7*r + 98*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**6*r - 640*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**10*r - 960*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**9*r - 752*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**8*r - 40*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**7*r - 98*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**6*r - 1740*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**11*r - 2172*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**10*r - 704*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**9*r + 584*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**8*r - 114*Kl_p*Ku_p*du**4*m**6*mu_l*mu_u*n**7*r - 680*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**12*r + 176*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**11*r + 1924*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**10*r + 1664*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**9*r + 576*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**13*r + 1008*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**12*r + 1564*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**11*r + 888*Kl_p*Ku_p*du**4*m**5*mu_l*mu_u*n**10*r - 50*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**14*r - 400*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**13*r - 128*Kl_p*Ku_p*du**4*m**4*mu_l*mu_u*n**12*r - 26*Kl_p*Ku_p*du**4*m**2*mu_l*mu_u*n**15*r - 88*Kl_p*Ku_p*du**4*m**3*mu_l*mu_u*n**14*r)/(1680*U_p**2*mu_l**2*mu_u**2*(m + n)**2*(m**2 + 4*m*n**3 + 6*m*n**2 + 4*m*n + n**4)**3)


def free_surface_J(**params):
    
    Ku_p, n, mu_u, m, U_p, n, du, F, r, Kl_p, m, mu_l = \
        (params["Ku_p"], params["n"], params["mu_u"], params["m"], params["U_p"], 
        params["n"], params["du"], params["F"], params["r"], params["Kl_p"], 
        params["m"], params["mu_l"])
    
    return -(700*Ku_p**2*du**4*m*mu_l*n**9 - 1260*Ku_p**2*du**4*mu_l*n**9 - 420*Ku_p**2*du**4*mu_l*n**10 - 2520*Ku_p**2*du**4*m*mu_l*n**7 - 1890*Ku_p**2*du**4*m*mu_l*n**8 - 945*Ku_p**2*du**4*mu_l*n**8 + 560*Ku_p**2*du**4*m*mu_l*n**10 + 378*Ku_p**2*du**4*mu_l*n**7*r + 1071*Ku_p**2*du**4*mu_l*n**8*r + 987*Ku_p**2*du**4*mu_l*n**9*r + 273*Ku_p**2*du**4*mu_l*n**10*r - 21*Ku_p**2*du**4*mu_l*n**11*r + 252*Ku_p**2*du**4*m**2*mu_l*n**4 + 96*Ku_p**2*du**4*m**3*mu_l*n**3 + 32*Ku_p**2*du**4*m**4*mu_l*n**2 + 756*Ku_p**2*du**4*m**2*mu_l*n**5 + 20*Ku_p**2*du**4*m**3*mu_l*n**4 - 32*Ku_p**2*du**4*m**5*mu_l*n**2 - 1232*Ku_p**2*du**4*m**2*mu_l*n**6 - 1220*Ku_p**2*du**4*m**3*mu_l*n**5 - 272*Ku_p**2*du**4*m**4*mu_l*n**4 - 96*Ku_p**2*du**4*m**5*mu_l*n**3 + 1288*Ku_p**2*du**4*m**2*mu_l*n**7 + 880*Ku_p**2*du**4*m**3*mu_l*n**6 + 464*Ku_p**2*du**4*m**4*mu_l*n**5 + 3059*Ku_p**2*du**4*m**2*mu_l*n**8 + 1232*Ku_p**2*du**4*m**3*mu_l*n**7 + 352*Ku_p**2*du**4*m**4*mu_l*n**6 + 784*Ku_p**2*du**4*m**2*mu_l*n**9 - 224*Ku_p**2*du**4*m**3*mu_l*n**8 - 140*Ku_p**2*du**4*m**2*mu_l*n**10 - 224*Ku_p**2*du**4*m**3*mu_l*n**9 + 672*Ku_p**2*du**4*m**2*mu_l*n**5*r + 1183*Ku_p**2*du**4*m**2*mu_l*n**6*r - 672*Ku_p**2*du**4*m**3*mu_l*n**5*r - 609*Ku_p**2*du**4*m**2*mu_l*n**7*r - 2065*Ku_p**2*du**4*m**3*mu_l*n**6*r - 1911*Ku_p**2*du**4*m**2*mu_l*n**8*r - 2058*Ku_p**2*du**4*m**3*mu_l*n**7*r - 763*Ku_p**2*du**4*m**2*mu_l*n**9*r - 672*Ku_p**2*du**4*m**3*mu_l*n**8*r + 42*Ku_p**2*du**4*m**2*mu_l*n**10*r + 3780*F**2*U_p**2*m**2*mu_l*mu_u**2*n**6 + 7560*F**2*U_p**2*m**3*mu_l*mu_u**2*n**5 + 3780*F**2*U_p**2*m**4*mu_l*mu_u**2*n**4 + 560*F**2*U_p**2*m**5*mu_l*mu_u**2*n**3 + 7560*F**2*U_p**2*m**2*mu_l*mu_u**2*n**7 + 12600*F**2*U_p**2*m**3*mu_l*mu_u**2*n**6 + 3360*F**2*U_p**2*m**4*mu_l*mu_u**2*n**5 + 6300*F**2*U_p**2*m**2*mu_l*mu_u**2*n**8 + 9240*F**2*U_p**2*m**3*mu_l*mu_u**2*n**7 + 1120*F**2*U_p**2*m**4*mu_l*mu_u**2*n**6 + 2520*F**2*U_p**2*m**2*mu_l*mu_u**2*n**9 + 3360*F**2*U_p**2*m**3*mu_l*mu_u**2*n**8 + 420*F**2*U_p**2*m**2*mu_l*mu_u**2*n**10 + 560*F**2*U_p**2*m**3*mu_l*mu_u**2*n**9 + 882*Ku_p**2*du**4*m*mu_l*n**6*r + 2289*Ku_p**2*du**4*m*mu_l*n**7*r + 1512*Ku_p**2*du**4*m*mu_l*n**8*r - 224*Ku_p**2*du**4*m*mu_l*n**9*r - 315*Ku_p**2*du**4*m*mu_l*n**10*r + 21*Ku_p**2*du**4*m*mu_l*n**11*r + 567*Kl_p*Ku_p*du**4*m**2*mu_u*n**7*r + 147*Kl_p*Ku_p*du**4*m**3*mu_u*n**6*r + 516*Kl_p*Ku_p*du**4*m**2*mu_u*n**8*r - 453*Kl_p*Ku_p*du**4*m**3*mu_u*n**7*r - 147*Kl_p*Ku_p*du**4*m**4*mu_u*n**6*r - 212*Kl_p*Ku_p*du**4*m**2*mu_u*n**9*r - 894*Kl_p*Ku_p*du**4*m**3*mu_u*n**8*r - 114*Kl_p*Ku_p*du**4*m**4*mu_u*n**7*r - 375*Kl_p*Ku_p*du**4*m**2*mu_u*n**10*r - 499*Kl_p*Ku_p*du**4*m**3*mu_u*n**9*r - 141*Kl_p*Ku_p*du**4*m**2*mu_u*n**11*r - 114*Kl_p*Ku_p*du**4*m**3*mu_u*n**10*r - 10*Kl_p*Ku_p*du**4*m**2*mu_u*n**12*r + 378*Kl_p*Ku_p*du**4*m*mu_u*n**8*r + 711*Kl_p*Ku_p*du**4*m*mu_u*n**9*r + 489*Kl_p*Ku_p*du**4*m*mu_u*n**10*r + 141*Kl_p*Ku_p*du**4*m*mu_u*n**11*r + 10*Kl_p*Ku_p*du**4*m*mu_u*n**12*r)/(1680*U_p**2*m**3*mu_l*mu_u**2*(n**3 + 3*n**2 + 3*n + m)**3)
