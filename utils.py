import numpy as np
import os
import matplotlib.pyplot as plt
import math

"""
utils.py by Sam Spinner. 
-------------------------
This file handles the majority of the graphing work needed to test the linear stability of superposed Plane Couette flow in 
the long wave perturbation approximation. 

The two functions to be used in a jupyter notebook are plot_vel_profile() and plot_stability_graph(). 
"""


def yihJ(**args):
    """ 
    Calculates J, equation (42) in Yih 1967.
    """ 

    # set up constants
    rho2, rho1, mu2, mu1, d2, d1, g, K, U0, verbose = \
        (args["rho2"], args["rho1"], args["mu2"], 
             args["mu1"], args["d2"], args["d1"], 
                 args["g"], args["dP"], args["U0"], args["verbose"])
    
    r = rho2 / rho1
    m = mu2 / mu1
    n = d2 / d1
    # these are all steady state quantities
    
    A2 = -(0.5 * K * U0) / mu2 * d1**2
    A1 = m * A2 

    a2 = (1 + (A2 * ((n ** 2) - m))) / (m + n)
    a1 = m * a2

    b = ((1 - A1 * (1 + n)) * n) / (m + n)

    # First alpha = 0 approx
    B1 = -((m + (3 * n ** 2) + (4 * n ** 3))/(2 * n ** 2 * (1 + n)))
    B2 = 2 * (m + (n ** 3))/(m * n) + (((n ** 2) * B1) / m)

    C2 = (m + n ** 3)/((m * n ** 2) * (1 + n))
    C1 = m * C2

    D2 = ((n ** 2) - m) / ((2 * m * n ** 2) * (1 + n))
    D1 = m * D2

    # this is for solution of upper boundary
    c_0p = (a2 - a1) / (B1 - B2)

    def h1(y):
        return (A1 * D1 * y ** 7) / 210 + \
            (a1 * D1 * y ** 6) / 60 + \
            ((a1 * C1 - (3 * c_0p * D1) - A1 * B1) * y ** 5) / 60 -\
            ((c_0p * C1 + A1) * (y ** 4) / 12)

    def h2(y):
        return (A2 * D2 * y ** 7) / 210 + \
            (a2 * D2 * y ** 6) / 60 + \
            ((a2 * C2 - 3 * c_0p * D2 - A2 * B2) * y ** 5) / 60 -\
            ((c_0p * C2 + A2) * (y ** 4) / 12)

    def h1p(y):
        return (A1 * D1 * y ** 6) / 30 + \
            (a1 * D1 * y ** 5) / 10 + \
            ((a1 * C1 - 3 * c_0p * D1 - A1 * B1) * y ** 4) / 12 -\
            ((c_0p * C1 + A1) * (y ** 3) / 3)

    def h2p(y):
        return (A2 * D2 * y ** 6) / 30 + \
            (a2 * D2 * y ** 5) / 10 + \
            ((a2 * C2 - 3 * c_0p * D2 - A2 * B2) * y ** 4) / 12 -\
            ((c_0p * C2 + A2) * (y ** 3) / 3)

    F2 = ((rho2 - rho1) / rho1) * (g * d1 / (U0 ** 2))

    h_1 = h1(1)
    h_1p = h1p(1)
    h_2 = h2(-n)
    h_2p = h2p(-n)
    
    if r == 1:
        H2 = r * h_2 
        J2 = r * h_2p 
    else:
        H2 = r * h_2 - ((n ** 3) / 6) * ((1 / ((c_0p * F2))) -(((r - 1) * (c_0p * B1 + a1))))
        J2 = r * h_2p + (((n ** 2) / 2)) * ((1 / ((c_0p * F2))) - (((r - 1) * (c_0p * B1 + a1))))
    
    
    J = (((1 / m) * (c_0p ** 2)) / (a1 - a2)) * \
    (m * (h_1p - 2 * h_1) - J2 - ((2 / n) * H2) + \
     (((m - n ** 2)/(2 * (1 + n))) * (h_1 - h_1p - (J2 / n) - (H2 / n ** 2))))
    
    if verbose:
        print('c_0p',c_0p)
        print(f'h_1: {h_1}, h_1p: {h_1p}, h_2: {h_2}, h2_p: {h_2p}')
        print(f'H2: {H2}, J2: {J2}')
        print('first term of J', (((1 / m) * (c_0p ** 2)) / (a1 - a2)))
        print('second term of J', (m * (h_1p - 2 * h_1) - J2 - ((2 / n) * H2) + (((m - n ** 2)/(2* (1 + n))) * (h_1 - h_1p - (J2 / n) - (H2 / n ** 2)))))
        print('J', J)
        print()
    return J

def calc_J_vectors(lines="d2s", variables="mu2s", **params):
    J_vectors = {}
    args = dict(params)

    for d2 in params[lines]:
        J_vector = []
        args[lines[:-1]] = d2

        for var in params[variables]:
            args[variables[:-1]] = var
            J_vector.append(yihJ(**args))

        J_vectors[d2] = J_vector
        
    return J_vectors


def generate_figure(J_vectors, x_axis, variable, path="images/", **params):

    if not os.path.exists(path):
        os.mkdir(path)
    
    if "figsize" not in params:
        params["figsize"] = (8,4)
    
    fig = plt.figure(figsize=params["figsize"])
    
    color_idx_red = 0
    color_idx_blue =0
    
    colors_red, colors_blue = (params["colors_red"], params["colors_blue"])
    
    for d2 in J_vectors:
        if d2 > 1:
            J = J_vectors[d2]
            plt.plot(params[x_axis], J, label=d2, color=colors_red[color_idx_red % len(colors_red)])
            color_idx_red += 1
        else:
            J = J_vectors[d2]
            plt.plot(params[x_axis], J, label=d2, color=colors_blue[color_idx_blue % len(colors_blue)])
            color_idx_blue += 1
    
    plt.axhline(y=0, color="black", linestyle="--")
    plt.legend()

    K = params["dP"]
    U0 = params["U0"]
    g = params["g"]
    m = round(params["mu2"] / params["mu1"], 3)
    r = round(params["rho2"] / params["rho1"], 3)
    
    if variable == "viscosity":
        plt.title("Stability of interface over changing {}, g,r,K,U0={},{},{},{}".format(variable, g, r, K, U0))
        plt.xlabel("m, μ2 / μ1")
        plt.xscale("log")
        
        if ("viscosity_lim" in params):
            plt.ylim(params["viscosity_lim"])
        
    if variable == "density":
        plt.title("Stability of interface over changing {}, g,m,K,U0={},{},{},{}".format(variable, g, m, K, U0))
        plt.xlabel("r, rho2 / rho1")
        plt.xscale("log")
        
        if ("density_lim" in params):
            plt.ylim(params["density_lim"])
        
    if variable == "differential pressure":
        plt.title("Stability of interface over changing {}, g,r,m,U0={},{},{},{}".format(variable, g, r, m, U0))
        plt.xlabel("K, - dP / dX")
        
        if ("pressure_lim" in params):
            plt.ylim(params["pressure_lim"])
        
    
    
    plt.ylabel("J, positive values represent instabilities")
    
    if "save" in params:
        filename = "{}_stability_".format(variable)
        
        if variable == "viscosity":
            filename += "r{}g{}K{}".format(r, g, K)
            
        if variable == "density":
            filename += "m{}g{}K{}".format(m, g, K)
            
        if variable == "differential pressure":
            filename += "m{}g{}r{}".format(m, g, r)
            
        if params["save"]:
            filename = params["save"]
        filename = filename.replace(".", ",")
        plt.savefig(path + filename)

    # When this function is called multiple times in one cell
    # there is a risk of overflowing ram. This usually happens at 
    # 20+ graphs opened simultaneously. In the case where one wants 
    # to generate many graphs at once but not necessarily see them within
    # cell of the jupyter notebook, uncomment the following line. 
        
    #plt.close(fig)

def plot_vel_profile(**params):
    """    
    Plots the dimensionless velocity profile of the superposed Plane Couette flow. 
    Velocity is made dimensionless by surface velocity U0. +y is "upper layer", 
    -y is "lower layer". 

    Equations (7) and (8) of Yih 1967. 
    """
    mu2, mu1, d2, d1, K, U0, N_POINTS = (params["mu2"], params["mu1"], params["d2"],
                                             params["d1"], params["dP"], params["U0"], params["N_POINTS"])
    
    # Calculate velocity profiles
    m = mu2 / mu1
    n = d2 / d1

    y1s = np.linspace(0, d1, N_POINTS) / d1
    y2s = np.linspace(-d2, 0, N_POINTS) / d1

    A2 = -(.5 * (K / mu2) * U0) * (d1 ** 2)
    A1 = m * A2

    a2 = (1 + A2 * (n**2 - m)) / (m + n)
    a1 = m * a2

    b = ((1 - A1*(1 + n)) * n) / (m + n)

    U1 = (A1 * (y1s ** 2) + a1 * y1s + b) / U0
    U2 = (A2 * (y2s ** 2) + a2 * y2s + b) / U0

    velocity_profile = np.concatenate((U2, U1))

    plt.figure(figsize=(8,4))
    
    # Visualize boundaries of layers and the interface
    plt.axhline(y=0, color="black", linestyle ="--") 
    plt.axhline(y= d1 / d1, color="black", linestyle ="--") 
    plt.axhline(y= -d2 / d1, color="black", linestyle ="--") 

    # generate plot, label axes, title
    plt.plot(velocity_profile, np.concatenate((y2s, y1s)))
    title="Velocity Profile n={}, m={}, K={}, U0={}".format(d2 / d1, mu2 / mu1, K, U0)
    plt.title(title)
    plt.ylabel("Flow height / d1, U1 > 0, U2 < 0")
    plt.xlabel("Velocity of flow")
    

def plot_stability_graph(**params):
    """
    Checks for the presence of "mu2s", "rho2s", and "dPs" within params. 
    For each one present, graphs the appropriate stability graph. 
    """
    
    # make viscosity change plot
    if "mu2s" in params.keys():
        J_vectors = calc_J_vectors(lines="d2s", variables="mu2s", **params)
        generate_figure(J_vectors, "mu2s", "viscosity", **params)
        
    # make density change plot
    if "rho2s" in params.keys():
        J_vectors = calc_J_vectors(lines="d2s", variables="rho2s", **params)
        generate_figure(J_vectors, "rho2s", "density", **params)
        
    # make K change plot
    if "dPs" in params.keys():
        J_vectors = calc_J_vectors(lines="d2s", variables="dPs", **params)
        generate_figure(J_vectors, "dPs", "differential pressure", **params)
