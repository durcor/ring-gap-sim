#!/usr/bin/env python3
import matplotlib.pyplot as plot
import csv
import subprocess
import sys

subprocess.run("make", check=True)
subprocess.run("./sim_ring_gaps " + ' '.join(sys.argv[1:]) + " > data.csv",
               check=True, shell=True)
data = csv.reader(open("data.csv"))
fig = plot.figure()
a = fig.add_axes([0, 0, 1, 1])
a.set_xlim(-300000000, 300000000)
a.set_ylim(-300000000, 300000000)
xm, ym = next(data)
plot.plot(float(xm), float(ym), "ro")
plot.plot(0, 0, "yo")
for x, y in data:
    plot.plot(float(x), float(y), "go", markersize=1)
plot.savefig("plot.png")
