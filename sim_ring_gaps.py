#!/usr/bin/env python
#
# Tyler Kaminski
# Adham Abdelwahab
# Vincent Andrade
# Professor Lukic
# PEP 351
#
# I pledge my honor that I have abided by the Stevens Honor System.
#
# Plot orbital positions to images using matplotlib.
import matplotlib.pyplot as plt
import numpy as np

# Constants for calculations
G = 6.67430e-11
AU = 1.495978707e11
# Mass of the central body (Saturn)
M = 5.6834e26
# Mass of Mimas
m = 3.7493e19
# total time of simulation (in seconds)
TIME = 100000000000000

particles = [[np.random.randint(-100, 100), np.random.randint(-100, 100)]
             for _ in range(1000)]


# General case for a single step in the RK# methods.
def get_ks(vx, vy, x, y, vxo, vyo, xo, yo, t):
    kx = (vx + vxo) * t
    ky = (vy + vyo) * t
    kvx = ((-G * M)/((x + xo)**2 + (y + yo)**2)) * (
        (x + xo)/np.sqrt((x + xo)**2 + (y + yo)**2)) * t
    kvy = ((-G * M)/((x + xo)**2 + (y + yo)**2)) * (
        (y + yo)/np.sqrt((x + xo)**2 + (y + yo)**2)) * t
    return kx, ky, kvx, kvy


# Simulate an orbit using a certain orbital approximation method.
def orbit(method, t, x, y, vx, vy):
    plt.plot(x, y, 'ro')
    for _ in range(TIME):
        # Simulate the euler, rk2, and rk4 methods based on input.
        if method == 'euler':
            r = np.sqrt(x**2 + y**2)

            ax = -G * M * (x/r**3)
            ay = -G * M * (y/r**3)

            x += vx * t + ax/2 * t**2
            y += vy * t + ay/2 * t**2

            vx += ax * t
            vy += ay * t
        elif method == 'rk2':
            k1x, k1y, k1vx, k1vy = get_ks(vx, vy, x, y,
                                          0, 0, 0, 0, t)
            k2x, k2y, k2vx, k2vy = get_ks(vx, vy, x, y,
                                          k1vx, k1vy, k1x, k1y, t)

            x += (k1x + k2x)/2
            y += (k1y + k2y)/2

            vx += (k1vx + k2vx)/2
            vy += (k1vy + k2vy)/2
        elif method == 'rk4':
            k1x, k1y, k1vx, k1vy = get_ks(vx, vy, x, y,
                                          0, 0, 0, 0, t)
            k2x, k2y, k2vx, k2vy = get_ks(vx, vy, x, y,
                                          k1vx/2, k1vy/2, k1x/2, k1y/2, t)
            k3x, k3y, k3vx, k3vy = get_ks(vx, vy, x, y,
                                          k2vx/2, k2vy/2, k2x/2, k2y/2, t)
            k4x, k4y, k4vx, k4vy = get_ks(vx, vy, x, y,
                                          k3vx, k3vy, k3x, k3y, t)

            x += (k1x + 2 * k2x + 2 * k3x + k4x)/6
            y += (k1y + 2 * k2y + 2 * k3y + k4y)/6

            vx += (k1vx + 2 * k2vx + 2 * k3vx + k4vx)/6
            vy += (k1vy + 2 * k2vy + 2 * k3vy + k4vy)/6
        # Plot each point and add its distance's absolute values to their
        # arrays so that eccentricity can be found.
        plt.plot(x, y, 'ro')


# Make the scales for the plots equal on both axes.
plt.gca().set_aspect('equal')

# Simulate each method using 29.8km/s, 42.1km/s, and 32.7km/s as inputs for
# the initial vy.
for method in ['euler', 'rk2', 'rk4']:
    # Start each orbit at Earth's location.
    t = TIME / 1e5
    vx = 0
    for p in particles:
        orbit(method, t, p[0], p[1], vx, vy=29.8e3)
    plt.savefig(method + ".png")
    plt.cla()
