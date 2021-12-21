#!/usr/bin/env python3
import matplotlib.pyplot as plot
import csv
import subprocess
import sys

if len(sys.argv) != 4:
    print("Usage:", sys.argv[0],
          "<num_particles> <mimas_multiplier> <num_time_steps>")
    exit(1)
subprocess.run("make", check=True)
subprocess.run("./sim_ring_gaps " + ' '.join(sys.argv[1:]) + " > data.csv",
               check=True, shell=True)
data = csv.reader(open("data.csv"))
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
for x, y in data:
    plot.plot(float(x), float(y), "go", markersize=1)
plot.legend(loc='upper right')
plot.savefig('img/mimas_mult=' + sys.argv[2] + '&p=' + sys.argv[1] + '&t=' +
             sys.argv[3] + '.png')
