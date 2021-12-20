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
m = 3.7493e19 * 10000
# Distance from Saturn to Mimas
d = 1.8552e8
# time step (length of one orbit of Mimas in seconds)
t = int(60 * 60 * 22.6 / 1000)
# total time of simulation (in mimas days)
TIME = 2000 * t


def orbit(interactive, num_particles, with_mimas, show_all):
    INNER_LIM = 100000000
    OUTER_LIM = 160000000

    # particles = np.random.randint(-OUTER_LIM, OUTER_LIM, (num_particles, 2))
    p = np.array([(float(i), 0.0) for i in range(
        INNER_LIM, OUTER_LIM, (OUTER_LIM - INNER_LIM) // num_particles)])
    v = np.array([(0.0, np.sqrt(G * M / x)) for x, _ in p])

    plot.gca().set_aspect('equal')

    # omega used for calculating position of mimas
    w = np.sqrt(G * M / d**3)
    xm, ym, rm = 0, 0, 0
    for n in range(0, TIME, t):
        if not show_all:
            plot.cla()
        if with_mimas:
            xm = d * np.cos(n * w)
            ym = d * np.sin(n * w)
        # Plot Saturn and Mimas
        if interactive or show_all or n == TIME - t:
            if with_mimas:
                plot.plot(xm, ym, 'ro')
            plot.plot(0, 0, "yo")
        for i in range(len(p)):
            rs = np.sqrt(p[i][0]**2 + p[i][1]**2)
            if with_mimas:
                rm = np.sqrt((p[i][0] - xm)**2 + (p[i][1] - ym)**2)
            if with_mimas:
                ax = -G * M * p[i][0] / rs**3 - G * m * (p[i][0] - xm) / rm**3
                ay = -G * M * p[i][1] / rs**3 - G * m * (p[i][1] - ym) / rm**3
            else:
                ax = -G * M * p[i][0] / rs**3
                ay = -G * M * p[i][1] / rs**3
            p[i][0] += v[i][0] * t + ax/2 * t**2
            p[i][1] += v[i][1] * t + ay/2 * t**2
            v[i][0] += ax * t
            v[i][1] += ay * t
            if interactive or show_all or n == TIME - t:
                plot.plot(p[i][0], p[i][1], "go", markersize=1)
        if interactive:
            plot.pause(0.00000001)
    if not interactive:
        plot.savefig('ring_gap.png')


orbit(interactive=False, num_particles=1000, with_mimas=True, show_all=True)
