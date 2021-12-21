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
plot.plot(float(xm), float(ym), "ro")
plot.plot(0, 0, "yo")
for x, y in data:
    plot.plot(float(x), float(y), "go", markersize=1)
plot.savefig('img/mimas_mult=' + sys.argv[2] + '&p=' + sys.argv[1] + '&t=' +
             sys.argv[3] + '.png')
