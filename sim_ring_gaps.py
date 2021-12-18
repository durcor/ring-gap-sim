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
m = 3.7493e19 * 1000
# Distance from Saturn to Mimas
d = 1.8552e8
# total time of simulation (in days)
TIME = 200
NUM_PARTICLES = 1000
LIM = 13000000000
method = 'euler'
t = 81360

# particles = np.random.randint(-LIM, LIM, (NUM_PARTICLES, 2))
p = np.array([[float(i), 0.0] for i in range(110000000, 130000000, 10000000)])


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


# TODO: Start Mimas at x=0
# Start all particles at x=0
v = np.array([[ 0, np.sqrt((G * M)/ x[0])] for x in p])
plot.scatter(0,0)
# Orbit simulation method to use. Avialable methods: 'euler', 'rk2', and 'rk4'
for n in range(0, TIME * t, t):
    xm = d * np.cos(n  * np.sqrt((G * M)/(d**3)))
    ym = d * np.sin(n  * np.sqrt((G * M)/(d**3)))
    plot.scatter(xm, ym)
    plot.pause(0.0000001)
    for i in range(len(p)):
        rs = np.sqrt(p[i][0]**2 + p[i][1]**2)
        rm = np.sqrt((p[i][0] - xm)**2 + (p[i][1] - ym)**2)

        ax = (-G * M * p[i][0]/rs**3) - ((G * m * (p[i][0] - xm)/rm**3))
        ay = (-G * M * p[i][1]/rs**3) - ((G * m * (p[i][1] - ym)/rm**3))

        p[i][0] += v[i][0] * t + ax/2 * t**2
        p[i][1] += v[i][0] * t + ay/2 * t**2
        v[i][0] += ax * t
        v[i][0] += ay * t
        ##plot.scatter(p[i][0], p[i][1])
        ##plot.pause(0.00000000000000000001)
plot.show()