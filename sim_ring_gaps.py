#!/usr/bin/env python3
import matplotlib.pyplot as plot
import csv
import subprocess
import sys
import numpy as np

if len(sys.argv) != 4:
    print("Usage:", sys.argv[0],
          "<num_particles> <mimas_multiplier> <num_time_steps>")
    exit(1)
subprocess.run(
    "clang -lm -fsanitize=safe-stack -Ofast -march=native sim_ring_gaps.c -o sim_ring_gaps",
    shell=True, check=True)
data_filename = 'data/mimas_mult=' + sys.argv[2] + '&p=' + sys.argv[1] + \
    '&t=' + sys.argv[3] + '.csv'
subprocess.run(
    "./sim_ring_gaps " + ' '.join(sys.argv[1:]) + " > '" + data_filename + "'",
    check=True, shell=True)
data = csv.reader(open(data_filename))
xm, ym = next(data)
plot.gca().set_aspect('equal')
plot.title("p=" + sys.argv[1] + ", m=" + sys.argv[2] +
           ", t=" + sys.argv[3])
plot.xlabel("Distance from Saturn (m)")
plot.ylabel("Distance from Saturn (m)")
plot.plot(float(xm), float(ym), "ro", label="Mimas")
plot.plot(0, 0, "yo", markersize=20, label="Saturn")
x, y = next(data)
plot.plot(float(x), float(y), "go", markersize=1, label="Ring Particle")
plot.legend(loc='upper right')

particles_in_division = 0.0
for x, y in data:
    if 1.1758e8 < np.sqrt(float(x)**2 + float(y)**2) < 1.2217e8:
        particles_in_division += 1
    plot.plot(float(x), float(y), "go", markersize=1)
print("Ratio of particles in division to total:",
      particles_in_division / int(sys.argv[1]))

plot.savefig('img/mimas_mult=' + sys.argv[2] + '&p=' + sys.argv[1] + '&t=' +
             sys.argv[3] + '.png')
