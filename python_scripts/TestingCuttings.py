import pyscan
import random
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from collections import deque
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


def plotInfTrapezoid(trapezoid, ax):

    low_x = max(trapezoid.left_x, ax.get_xlim()[0])
    high_x = min(trapezoid.right_x, ax.get_xlim()[1])

    def boundU(t):
        if t == float('inf'):
            return ax.get_ylim()[1]
        return t
    def boundL(l):
        if l == -float('inf'):
            return ax.get_ylim()[0]
        return l

    tRy = min(trapezoid.top_a * high_x + trapezoid.top_b, ax.get_ylim()[1])
    tLy = min(trapezoid.top_a * low_x + trapezoid.top_b, ax.get_ylim()[1])
    bRy = max(trapezoid.bottom_a * high_x + trapezoid.bottom_b, ax.get_ylim()[0])
    bLy = max(trapezoid.bottom_a * low_x + trapezoid.bottom_b, ax.get_ylim()[0])
    verts = [(low_x, tLy),
            (high_x, tRy),
            (high_x, bRy),
            (low_x, bLy),
            (low_x, tLy)]
    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, lw=2, facecolor='orange', alpha=.2)
    ax.add_patch(patch)
    #print(trapezoid)



def plotLine(line, x_low, x_high, ax):
    ax.plot([x_low, x_high], [line[0] * x_low + line[1], line[0] * x_high + line[1]], color='k')


def testCrossing(trap, line):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plotTrapezoid(trap, ax)

    ys = [trap.top_a * trap.left_x + trap.top_b,
      trap.top_a * trap.right_x + trap.top_b,
      trap.bottom_a * trap.right_x + trap.bottom_b,
      trap.bottom_a * trap.left_x + trap.bottom_b,
      trap.top_a * trap.left_x + trap.top_b]
    ax.set_ylim([min(*ys), max(*ys)])
    ax.set_xlim([trap.left_x - 1, trap.right_x + 1])

    plotLine(line, trap.left_x, trap.right_x, ax)
    traps = pyscan.splitCell(trap, [], line)
    for t in traps:
        plotTrapezoid(t, ax)
    plt.show()

def testCrossesTopBottom():
    trap = pyscan.Trapezoid(-0.214544,-0.167938, -0.789927,  0.789774, -0.961276,  -0.324445,)
    line = (10, 4.2)
    testCrossing(trap, line)

def testCrossesTopBottomDegPos():

    trap = pyscan.Trapezoid(1,0, pyscan.xLoc(1, 0, .2, .1), .2, .1, 1)
    line = (2.5, -.9)
    testCrossing(trap, line)

def testCrossesTopBottomDegNeg():
    trap = pyscan.Trapezoid(1,0, pyscan.xLoc(1, 0, .2, .1), .2, .1, 1)
    line = (-2.5, 1.8)
    testCrossing(trap, line)

def testCrossesLeftTop():
    trap = pyscan.Trapezoid(-0.214544,-0.167938, -0.789927,  0.789774, -0.961276,  -0.324445,)
    line = (.8, .5)
    testCrossing(trap, line)


def testCrossesRightTop():
    trap = pyscan.Trapezoid(-0.214544,-0.167938, -0.789927,  0.789774, -0.961276,  -0.324445,)
    line = (-.8, -.5)
    testCrossing(trap, line)

def testCrossesLeftRight():
    trap = pyscan.Trapezoid(-0.214544,-0.167938, -0.789927,  0.789774, -0.961276,  -0.324445)
    line = (.8, -.7)
    testCrossing(trap, line)

def testCrossesLeftBottom():
    trap = pyscan.Trapezoid(-0.214544,-0.167938, -0.789927,  0.789774, -0.961276,  -0.324445)
    line = (-.9, -1.9)
    testCrossing(trap, line)

def testCrossesRightBottom():
    trap = pyscan.Trapezoid(-0.214544,-0.167938, -0.789927,  0.789774, -0.961276,  -0.324445)
    line = (2, -.4)
    testCrossing(trap, line)

def testInfTrap(trap, line, low_x, high_x):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim([low_x, high_x])
    ax.set_ylim([low_x, high_x])
    plotInfTrapezoid(trap, ax)
    traps = pyscan.splitCell(trap, [], line)
    for t in traps:
        plotInfTrapezoid(t, ax)
    plotLine(line, low_x, high_x, ax)
    plt.show()
    #testCrossing(trap, line)

    plt.show()

def testInfR_TR():
    trap = pyscan.Trapezoid(0.79867, -0.384848, -0.576201, 0.345627, -0.645891, float("inf"))
    line = (.5, -.4)
    testInfTrap(trap, line, -2, 2)

def testInfR_LR():
    trap = pyscan.Trapezoid(0.79867, -0.1, -0.576201, 0.345627, -0.645891, float("inf"))
    line = (.5, -.4)
    testInfTrap(trap, line, -2, 2)


def testInfR_BR():
    trap = pyscan.Trapezoid(0.79867, -0.1, pyscan.xLoc(.79867, -0.1, .345627, -.345891), 0.345627, -0.345891, float("inf"))
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #plotInfTrapezoid(trap, ax, -2, 2)
    #plt.show()
    line = (.5, -.4)
    testInfTrap(trap, line, -2, 2)

def testInfR_BL():
    trap = pyscan.Trapezoid(0.345627, -0.345891, -float('inf'), 0.79867, -0.1, pyscan.xLoc(.79867, -0.1, .345627, -.345891))
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #plotInfTrapezoid(trap, ax, -2, 2)
    #plt.show()
    line = (.5, -.4)
    testInfTrap(trap, line, -2, 2)

def testInfR_TL():
    trap = pyscan.Trapezoid(0.345627, -0.345891, -float('inf'), 0.79867, -0.1, pyscan.xLoc(.79867, -0.1, .345627, -.345891))
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #plotInfTrapezoid(trap, ax, -2, 2)
    #plt.show()
    line = (.5, -.2)
    testInfTrap(trap, line, -2, 2)


def testInf():
    trap = pyscan.Trapezoid(0, float('inf'), -float('inf'), 0, -float('inf'), float('inf'))
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #plotInfTrapezoid(trap, ax, -2, 2)
    #plt.show()
    line = (.5, -.2)
    testInfTrap(trap, line, -2, 2)

def testInfLower():
    trap = pyscan.Trapezoid(0.5, -0.2, -float('inf'), 0, -float('inf'), float('inf'))
    line = (.4, -.2)
    testInfTrap(trap, line, -2, 2)

def testInfHigher():
    trap = pyscan.Trapezoid(0, float('inf'), -float('inf'), 0.5, -0.2, float('inf'))
    line = (.4, -.2)
    testInfTrap(trap, line, -2, 2)


def testInfHigherLeft():
    trap = pyscan.Trapezoid(0, float('inf'), 0, 0.2, -0.2, float('inf'))
    line = (.5, -.3)
    testInfTrap(trap, line, -2, 2)

    line = (.5, -.1)
    testInfTrap(trap, line, -2, 2)

    line = (.1, -.1)
    testInfTrap(trap, line, -2, 2)


def testInfLowerLeft():
    trap = pyscan.Trapezoid(0.5, -0.2, 0, 0, float("-inf"), float('inf'))
    line = (.2, -.4)
    testInfTrap(trap, line, -2, 2)

    line = (.7, -.5)
    testInfTrap(trap, line, -2, 2)
    line = (.2, -.13)
    testInfTrap(trap, line, -2, 2)


def testInfHigherRight():
    trap = pyscan.Trapezoid(0, float('inf'), float('-inf'), 0.4, -0.2, 0)
    line = (.5, -.3)
    testInfTrap(trap, line, -2, 2)

    line = (.5, -.1)
    testInfTrap(trap, line, -2, 2)

    line = (.1, -.1)
    testInfTrap(trap, line, -2, 2)


def testInfLowerRight():
    trap = pyscan.Trapezoid(0.5, -0.2, float('-inf'), 0, float('-inf'), 0)
    line = (.2, -.4)
    testInfTrap(trap, line, -2, 2)

    line = (.7, -.5)
    testInfTrap(trap, line, -2, 2)
    line = (.2, -.13)
    testInfTrap(trap, line, -2, 2)

def plotTrap(trap, low_x, high_x):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plotInfTrapezoid(trap, ax, low_x, high_x)
    plt.show()


"""
testCrossesTopBottomDegPos()
testCrossesTopBottomDegNeg()

testCrossesLeftTop()
testCrossesRightTop()
testCrossesLeftBottom()
testCrossesRightBottom()

testCrossesLeftRight()
testInfR_BR()
testInfR_BL()

testInfR_TR()
testInfR_TL()


testInf()
testInfLower()
testInfHigher()
testInfLowerLeft()
testInfHigherLeft()
testInfLowerRight()
testInfHigherRight()

"""

def plotPartitioning(points, cells):
    x, y = zip(*points)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, '.')
    for cell in cells:
        plotInfTrapezoid(cell[0], ax)
    plt.show()

def halfPartitioning(points, r, sample_size):
    cells = deque()
    cells.append(points)
    final_cells = []
    while len(cells) + len(final_cells) < sample_size:

        cell_points = cells.popleft()
        print len(cell_points)
        if len(cell_points) <= len(points) / sample_size:
            final_cells.append(cell_points)
            continue
        trap_lists = pyscan.createPartition(cell_points, r)
        #plotPartitioning(cell_points, trap_lists)

        left_overs = set(cell_points)
        for inner_cell_points in trap_lists:
            left_overs.difference_update(inner_cell_points[2])
            cells.append(inner_cell_points[2])
        print list(left_overs)
        cells.append(list(left_overs))
    print list(cells) + final_cells
    return list(cells) + final_cells

def partitioningSample(pointCells):
    sample = []
    for cell in pointCells:
        chosen_int = random.randint(0, len(cell) - 1)
        if (len(cell) == 3):
            print cell
        sample.append((cell[chosen_int][0], cell[chosen_int][1], float(len(cell))))
    return sample


def plotWeightedSample(sample):
    x, y, c = zip(*sample)
    plt.scatter(x, y, c=c, cmap='gray')
    plt.show()
    #plt.gray()

if __name__ == "__main__":
    line_count = 200
    sample_size = 20
    r = 4
    count = 100
    low_x = -2
    high_x = 2


    linelist = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(line_count)]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    weightedSample = partitioningSample(halfPartitioning(linelist, r, sample_size))
    print weightedSample
    plotWeightedSample(weightedSample)


    if True:
        x, y = zip(*linelist)
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        ax.plot(x, y, '.')
        cells = pyscan.createPartition(linelist, r)
        
        for cell in cells:
           plotInfTrapezoid(cell[0], ax)
        test_set = pyscan.testSet(linelist, r)
        #cells = pyscan.randomCuttings(test_set, r)
        #for line in test_set:
        #    plotLine(line, low_x, high_x, ax)

        ax.set_xlim([low_x, high_x])
        ax.set_ylim([low_x, high_x])
    else:

        cells = pyscan.randomCuttings(linelist, r)
        print len(cells)
        for line in linelist:
            plotLine(line, low_x, high_x, ax)
        for cell in cells:
            plotInfTrapezoid(cell, ax)
    #plt.savefig('destination_path.svg', format='svg', dpi=10000)
    plt.show()
    """
    num_cells = 0
    for i in xrange(count):
        linelist = [(2 * random.random() - 1, 2 * random.random() - 1) for i in range(line_count)]
        cells = pyscan.randomCuttings(linelist, r)
        num_cells += len(cells)
    print (num_cells / float(count)) / (r * r)
    """

