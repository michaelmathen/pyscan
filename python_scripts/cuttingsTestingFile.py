import pyscan
import unittest
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

def plotTrapezoid(trapezoid, ax, color='orange'):
    print(trapezoid)
    verts = [(trapezoid.left_x, trapezoid.top_a * trapezoid.left_x + trapezoid.top_b),
             (trapezoid.right_x, trapezoid.top_a * trapezoid.right_x + trapezoid.top_b),
             (trapezoid.right_x, trapezoid.bottom_a * trapezoid.right_x + trapezoid.bottom_b),
             (trapezoid.left_x, trapezoid.bottom_a * trapezoid.left_x + trapezoid.bottom_b),
             (trapezoid.left_x, trapezoid.top_a * trapezoid.left_x + trapezoid.top_b)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, lw=2, facecolor=color, alpha=.2)
    ax.add_patch(patch)

def plotLine(line, x_low, x_high, ax):
    ax.plot([x_low, x_high], [line[0] * x_low + line[1], line[0] * x_high + line[1]])

class TestSplitCell(unittest.TestCase):

    def testBottomRightSplit(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        line = (-0.314999, -0.259759)
        low_x = -2
        high_x = 2
        trap1 = pyscan.Trapezoid(0.0201098, 0.264799, -0.412858,  -0.721641, -0.196454, 0.621845)

        trap2 = pyscan.Trapezoid(0.0201098, 0.264799, -0.412858, -0.721641, -0.196454, 0.155677)
        trap3 = pyscan.Trapezoid(0.0201098, 0.264799, 0.155677, -0.314999, -0.259759, 0.621845)
        trap4 = pyscan.Trapezoid(-0.314999, -0.259759, 0.155677, -0.721641, -0.196454, 0.621845)
        cells = pyscan.splitCell(trap1, [], (-0.314999, -0.259759))
        self.assertEqual(len(cells), 3)
        self.assertEqual(trap2, cells[0])
        self.assertEqual(trap3, cells[1])
        self.assertEqualt(trap4, cells[2])
        plotLine(line, low_x, high_x, ax)
        plotTrapezoid(trap1, ax)

        for trap in pyscan.splitCell(trap1, [], (-0.314999, -0.259759)):
            plotTrapezoid(trap, ax, color='blue')

        plt.show()


def plotTrap1(line):

if __name__ == '__main__':

    trap1 = pyscan.Trapezoid(0.0201098, 0.264799, -0.412858,  -0.721641, -0.196454, 0.621845)
    line = (-0.314999, -0.259759)
    unittest.main()

