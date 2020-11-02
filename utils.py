import numpy as np
import os
import matplotlib.pyplot as plt
import math
from itertools import cycle

"""
utils.py by Sam Spinner. 
-------------------------
This file handles the majority of the graphing work needed to test the linear stability of superposed Plane Couette flow in 
the long wave perturbation approximation. 

The two functions to be used in a jupyter notebook are plot_vel_profile() and plot_stability_graph(). 
"""


def yihJ(**params):
    """ 
    Calculates J, equation (42) in Yih 1967.
    """ 
    # set up constants
    rho2, rho1, mu2, mu1, d2, d1, g, K, U0, verbose = \
        (params["rho2"], params["rho1"], params["mu2"], 
             params["mu1"], params["d2"], params["d1"], 
                 params["g"], params["dP"], params["U0"], params["verbose"])
    
    r = rho2 / rho1
    m = mu2 / mu1
    n = d2 / d1

    args = (r, m, n, g, K, U0) 

    cache = False
    if "cache" in params:
        cache = True
   
    if cache and args in params["cache"]:
        cache = False
        return params["cache"]
    
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
    
    H2 = r * h_2 - ((n ** 3) / 6) * ((1 / ((c_0p))) * F2) -(((r - 1) * (c_0p * B1 + a1)))
    J2 = r * h_2p + (((n ** 2) / 2)) * ((1 / ((c_0p))) * F2) - (((r - 1) * (c_0p * B1 + a1)))
    
    
    J = (((1 / m) * (c_0p ** 2)) / (a1 - a2)) * \
    (m * (h_1p - 2 * h_1) - J2 - ((2 / n) * H2) + \
     (((m - n ** 2)/(2 * (1 + n))) * (h_1 - h_1p - (J2 / n) - (H2 / n ** 2))))
    math
    if verbose:
        print('c_0p',c_0p)
        print(f'h_1: {h_1}, h_1p: {h_1p}, h_2: {h_2}, h2_p: {h_2p}')
        print(f'H2: {H2}, J2: {J2}')
        print('first term of J', (((1 / m) * (c_0p ** 2)) / (a1 - a2)))
        print('second term of J', (m * (h_1p - 2 * h_1) - J2 - ((2 / n) * H2) + (((m - n ** 2)/(2* (1 + n))) * (h_1 - h_1p - (J2 / n) - (H2 / n ** 2)))))
        print('J', J)
        print()

    if cache:
        params["cache"][args] = J

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


def generate_figure(J_vectors, x_axis, variable, ax=None, path="images/", **params):

    if not os.path.exists(path):
        os.mkdir(path)
    
    if "stability_figsize" not in params:
        params["stability_figsize"] = (8,4)
    
    ax = ax or plt.gca()

    # remove the top, bottom, and left axes
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    color_idx_red = 0
    color_idx_blue =0
    
    colors_red, colors_blue = (params["colors_red"], params["colors_blue"])

    
    if "label_size" not in params:
        params["label_size"] = 14
    if "title_size" not in params:
        params["title_size"] = 14

    label_size = params["label_size"]
    title_size = params["title_size"]
    
    for d2 in J_vectors:
        if d2 > 1:
            J = J_vectors[d2]
            ax.plot(params[x_axis], J, label=d2, color=colors_red[color_idx_red % len(colors_red)])
            color_idx_red += 1
        elif d2 < 1:
            J = J_vectors[d2]
            ax.plot(params[x_axis], J, label=d2, color=colors_blue[color_idx_blue % len(colors_blue)])
            color_idx_blue += 1
        # top and bottom layer are same height
        else:
            J = J_vectors[d2]
            ax.plot(params[x_axis], J, label=d2, color=params["color_one"])
            
    
    ax.axhline(y=0, color="gray", linestyle="--")

    # we don't want to automatically include a legend
    # if user has made legend settings

    K = params["dP"]
    U0 = params["U0"]
    g = params["g"]
    m = round(params["mu2"] / params["mu1"], 3)
    r = round(params["rho2"] / params["rho1"], 3)

    include_legend = False
    
    if variable == "viscosity":
        if "visc_title" in params and params["visc_title"] is not None:
            title = params["visc_title"]
        else: 
            title = "Stability of interface over changing {}, g,r,k,u0={},{},{},{}".format(variable,
                    g, r, K, U0)
        if "x_label" in params and params["x_label"] is not None:
            xlabel = params["x_label"]
        else:
            xlabel = "$m \ \ (\mu_2 / \mu_1)$"
        ax.set_xscale("log")
        
        if ("viscosity_lim" in params):
            ax.set_ylim(params["viscosity_lim"])

        if "legends" in params and params["legends"][0]:
            include_legend = True
        
    if variable == "density":

        if "dens_title" in params:
            title = params["dens_title"]
        else: 
            title = "Stability of interface over changing {}, g,m,K,U0={},{},{},{}".format(variable, g, m, K, U0)

        # unicode character for lowercase rho
        xlabel = "$r \ \  (\u03C1_2 / \u03C1_1)$"
        ax.set_xscale("log")
        
        if ("density_lim" in params):
            ax.set_ylim(params["density_lim"])
        
        if "legends" in params and params["legends"][1]:
            include_legend = True

    if variable == "differential pressure":

        if "pressure_title" in params:
            title = params["pressure_title"]
        else: 
            title = "Stability of interface over changing {}, g,r,m,U0={},{},{},{}".format(variable, g, r, m, U0)
        xlabel = "$K \ \ (-\partial p / \partial X)$"
        
        if ("pressure_lim" in params):
            ax.set_ylim(params["pressure_lim"])

        if "legends" in params and params["legends"][2]:
            include_legend = True
    
    if "legends" not in params or include_legend:
        ax.legend(title="$ n \ \ (d_2 / d_1)$")
    
    ax.set_title(title, fontsize=title_size) 
    ax.set_xlabel(xlabel, fontsize=label_size) 
    ax.set_ylabel("J", fontsize=label_size)
    

def plot_vel_profile(ax, path="images/", **params):
    """    
    Plots the dimensionless velocity profile of the superposed Plane Couette flow. 
    Velocity is made dimensionless by surface velocity U0. +y is "upper layer", 
    -y is "lower layer". 

    Equations (7) and (8) of Yih 1967. 
    """
    mu2s, mu1, d2, d1, K, U0, N_POINTS = (params["mu2s"], params["mu1"], params["d2"],
                                             params["d1"], params["dP"], params["U0"], params["N_POINTS"])
    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)

    colorcycler = None
    if "vel_colors" in params and params["vel_colors"] is not None:
        colorcycler = cycle(params["vel_colors"]) 

    for mu2 in mu2s:
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

        # generate plot
        if colorcycler:
            ax.plot(velocity_profile, np.concatenate((y2s, y1s)), label=f"{mu2 / mu1}", 
                    linestyle=next(linecycler), color=next(colorcycler))
        else:
            ax.plot(velocity_profile, np.concatenate((y2s, y1s)), label=f"{mu2 / mu1}", 
                    linestyle=next(linecycler))
            

    if "legend" in params and params["legend"]:
        ax.legend(title="$m$")

    # remove the top, bottom, and left axis lines
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # Visualize boundaries of layers and the interface
    ax.axhline(y=0, color="gray", linestyle ="--", linewidth=1) 
    ax.axhline(y= d1 / d1, color="black", linewidth=1) 
    ax.axhline(y= -d2 / d1, color="black", linewidth=1) 

    # Create title and axes labels
    if "vel_title" in params and params["vel_title"] is not None:
        title = params["vel_title"]
    else:
        title="Velocity Profile n={}, m={}, K={}, U0={}".format(d2 / d1, mu2 / mu1, K, U0)

    if "label_size" not in params:
        params["label_size"] = 14
    if "title_size" not in params:
        params["title_size"] = 14
    
    # Title and label axes
    ax.set_title(title, fontsize=params["title_size"])

    if "y_label" in params and params["y_label"] is not None:
        ax.set_ylabel(params["y_label"], fontsize=params["label_size"])
    if "x_label" in params and params["x_label"] is not None:
        ax.set_xlabel(params["x_label"], fontsize=params["label_size"])


def plot_stability_graph(ax, arg="viscosity", **params):
    """
    Checks for the presence of "mu2s", "rho2s", and "dPs" within params. 
    For each one present, graphs the appropriate stability graph. 
    """
    if arg not in ["viscosity", "density", "velocity"]:
        raise ValueError("Incorrect arg for plot_stability_graph: {}".format(arg))
    
    # make viscosity change plot
    if arg == "viscosity":
        J_vectors = calc_J_vectors(lines="d2s", variables="mu2s", **params)
        generate_figure(J_vectors, "mu2s", "viscosity", ax, **params)
        
    # make density change plot
    if arg == "density":
        J_vectors = calc_J_vectors(lines="d2s", variables="rho2s", **params)
        generate_figure(J_vectors, "rho2s", "density", ax, **params)
        
    # make K change plot
    if arg == "velocity":
        J_vectors = calc_J_vectors(lines="d2s", variables="dPs", **params)
        generate_figure(J_vectors, "dPs", "differential pressure", ax, **params)
