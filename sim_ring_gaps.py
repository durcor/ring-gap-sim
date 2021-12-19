#!/usr/bin/env python3
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
import matplotlib.pyplot as plot
import numpy as np

# Constants for calculations
G = 6.67430e-11
AU = 1.495978707e11
# Mass of the central body (Saturn)
M = 5.6834e26
# Mass of Mimas
m = 3.7493e19
# Distance from Saturn to Mimas
d = 1.8552e8
method = 'euler'
# time step (length of one orbit of Mimas in seconds)
t = int(60 * 60 * 22.6 / 1000)
# total time of simulation (in mimas days)
TIME = 2000 * t

INNER_LIM = int(d - 100000000)
# INNER_LIM = 100000000
OUTER_LIM = int(d + 100000000)
# OUTER_LIM = 160000000
NUM_PARTICLES = 10

# particles = np.random.randint(-LIM, LIM, (NUM_PARTICLES, 2))
p = np.array([(float(i), 0.0) for i in range(
    INNER_LIM, OUTER_LIM, (OUTER_LIM - INNER_LIM) // NUM_PARTICLES)])
v = np.array([(0.0, np.sqrt(G * M / x)) for x, _ in p])


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
    # Simulate the euler, rk2, and rk4 methods based on input.
    if method == 'euler':
        rs = np.sqrt(x**2 + y**2)
        rm = np.sqrt((x - xm)**2 + (y - ym)**2)

        ax = -G * M * (x/rs**3) - G * m * (x - xm)/rm**3
        ay = -G * M * (y/rs**3) - G * m * (y - ym)/rm**3

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
    plot.scatter(x, y)


plot.gca().set_aspect('equal')

# omega used for calculating position of mimas
w = np.sqrt(G * M / d**3)
for n in range(0, TIME, t):
    plot.cla()
    xm = d * np.cos(n * w)
    ym = d * np.sin(n * w)
    plot.plot(xm, ym, 'ro')
    # Plot Saturn
    plot.plot(0, 0, "yo")
    for i in range(len(p)):
        rs = np.sqrt(p[i][0]**2 + p[i][1]**2)
        rm = np.sqrt((p[i][0] - xm)**2 + (p[i][1] - ym)**2)

        ax = -G * M * p[i][0] / rs**3 - G * m * (p[i][0] - xm) / rm**3
        ay = -G * M * p[i][1] / rs**3 - G * m * (p[i][1] - ym) / rm**3

        p[i][0] += (v[i][0] * t) + (ax/2 * t**2)
        p[i][1] += (v[i][1] * t) + (ay/2 * t**2)
        v[i][0] += ax * t
        v[i][1] += ay * t
        plot.plot(p[i][0], p[i][1], "go")
    # plot.pause(0.00000001)
plot.savefig('ring_gap.png')
