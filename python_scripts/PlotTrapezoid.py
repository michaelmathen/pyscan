import pyscan
import random
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np


def plotTrapezoid(trapezoid, ax):
    print(trapezoid)
    verts = [(trapezoid.left_x, trapezoid.top_a * trapezoid.left_x + trapezoid.top_b),
             (trapezoid.right_x, trapezoid.top_a * trapezoid.right_x + trapezoid.top_b),
             (trapezoid.right_x, trapezoid.bottom_a * trapezoid.right_x + trapezoid.bottom_b),
             (trapezoid.left_x, trapezoid.bottom_a * trapezoid.left_x + trapezoid.bottom_b),
             (trapezoid.left_x, trapezoid.top_a * trapezoid.left_x + trapezoid.top_b)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, lw=2, facecolor='orange', alpha=.2)
    ax.add_patch(patch)

def plotLine(line, x_low, x_high, ax):
    ax.plot([x_low, x_high], [line[0] * x_low + line[1], line[0] * x_high + line[1]])




low_x = -2
high_x = 2

fig = plt.figure()
ax = fig.add_subplot(111)
line = (-0.314999, -0.259759)
cell = pyscan.Trapezoid(-0.721641, -0.196454, -0.412858, 0.0201098, 0.264799, 0.621845)

plotTrapezoid(cell, ax)

print("left = " + str(cell.crossesLeftSide(line[0], line[1])))
print("right = " + str(cell.crossesRightSide(line[0], line[1])))
print("top = " + str(cell.crossesTop(line[0], line[1])))
print("bottom = " + str(cell.crossesBottom(line[0], line[1])))
plotLine(line, low_x, high_x, ax)

plt.show()