import pyscan
import argparse
import random
import numpy.random as npr
import csv
import bisect
import time
import math
import matplotlib.pyplot as plt
import eps_scan
import sys

def normalize(a, b):
    mag = math.sqrt(a*a + b*b)
    return (a / mag, b / mag)

"""
Given a crime keyword find a region that is interesting and display that region. 
Also explain how interesting that region is.
"""

def randomLinFunc():
    a = 0
    b = 0
    while not ((a < 0 < b) or (b < 0 < a)):
        a = npr.randn()
        b = npr.randn()
    mag = math.sqrt(a*a + b*b)
    print a, b
    return (a / mag, b / mag)

def permSampleL(collection, s, n):
    s_a = random.sample(collection, s)
    s_b = random.sample(collection, s)
    N = random.sample(collection, n)
    return (N, s_a, s_b)

def mkCDF(statistic):
    values = sorted(statistic)
    return lambda x: float(bisect.bisect_left(values, x) - 1) / len(values)


def actualSubgridValue(grid, subgrid, a, b):
    gridVal = 0
    for i in xrange(subgrid.upCol(), subgrid.lowCol() + 1):
        for j in xrange(subgrid.upRow(), subgrid.lowRow() + 1):
            gridVal += a * grid.redWeight(j, i) + b * grid.blueWeight(j, i)
    return gridVal

def toEpsScanPoints(points):
    return [eps_scan.Point(x, y, val) for x, y, val in points]

def toPyScanPoints(points):
    return [pyscan.Point(p.x, p.y, p.anomaly, not p.anomaly) for p in points]

def toTuples(points):
    return [(p.getX(), p.getY(), p.getRedWeight()) for p in points]

def randomRect(points, r):
    def bBox(*bb):
        return (min(*map(lambda pt: pt.x, bb)),
                max(*map(lambda pt: pt.x, bb)),
                min(*map(lambda pt: pt.y, bb)),
                max(*map(lambda pt: pt.y, bb)))
    while True:
        [pt1, pt2] = random.sample(points, 2)

        if random.randint(0, 1) == 0:
            (minX, maxX, minY, maxY) = bBox(pt1, pt2)
            fx = lambda pt: pt.x
            fy = lambda pt: pt.y
        else:
            fx = lambda pt: pt.y
            fy = lambda pt: pt.x
            (minY, maxY, minX, maxX) = bBox(pt1, pt2)

        if random.randint(0, 1) == 0:
            sf = lambda pt: fx(pt)
        else:
            sf = lambda pt: -fx(pt)
            minX, maxX = maxX, minX

        subPoints = [pt for pt in points if minY <= fy(pt) and fy(pt) <= maxY]
        # Make a pass through the data and find the lower and upper bound.
        inB = [pt for pt in subPoints if minX <= fx(pt) and fx(pt) <= maxX]
        if random.randint(0, 1):  # random.randint(0, 1) == 0:
            ab = [pt for pt in subPoints if maxX < fx(pt)]
            ab.sort(key=lambda p: -sf(p))
        else:
            ab = [pt for pt in subPoints if fx(pt) < minX]
            ab.sort(key=lambda p: sf(p))
        while len(ab) > 0:
            if len(inB) > int(r * len(points) + .5):
                break
            if len(inB) == int(r * len(points) + .5):
                return eps_scan.Rectangle(*bBox(*inB))
            el = ab.pop()
            inB.append(el)
            while len(ab) > 0 and sf(el) == sf(ab[-1]):
                el = ab.pop()
                inB.append(el)


def maxSymDiff(A, B, points, grids):
    B = grids.toRectangle(B)
    A_Cap_B = 0
    A_Cup_B = 0
    print A
    print B
    for pt in points:
        if A.inside(eps_scan.Point(pt.getX(), pt.getY(), pt.getRedWeight())) or B.contains(pt):
            A_Cup_B += 1
        if A.inside(eps_scan.Point(pt.getX(), pt.getY(), pt.getRedWeight())) and B.contains(pt):
            A_Cap_B += 1
    return 1 - (A_Cap_B) / float(A_Cup_B)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--n', type=int,
                        default=20,
                        help='The size of the grid')
    parser.add_argument("--a", type=float,
                        default=.1,
                        help='red weight')
    parser.add_argument("--b", type=float,
                        default=-.1,
                        help="blue weight")
    parser.add_argument('--repr', type=int,
                        default=12,
                        help='The max grid approximation.')
    parser.add_argument('--repeats', type=int,
                        default=20,
                        help='Number of time experiment is repeated.')
    parser.add_argument('--disable_slow',
                        action='store_true',
                        help="Disable the slow algorithm")
    parser.add_argument('--prefix',
                        required=True,
                        help="prefix for the timing and error outputs")
    parser.add_argument('--r', type=float,
                        default=.2,
                        help="Max region size")
    parser.add_argument('--red_bias', type=float,
                        default=.1,
                        help="Increased bias of red points in max region")

    parser.add_argument('--max_time',
                        type=float,
                        default=5,
                        help='The maximum time we allow')

    args = parser.parse_args()

    (a, b) = normalize(args.a, args.b)
    args.a = a
    args.b = b
    csv_points = []
    fieldnames = ["time", "maxDiff", "alg", "r"]
    time_file = open(args.prefix + "_etime.csv", 'w')
    timingWriter = csv.DictWriter(time_file, fieldnames=fieldnames)
    timingWriter.writeheader()
    with open("police_inct.csv", "r") as f:
        reader = csv.DictReader(f)
        points = []
        numAnomalies = 0
        for row in reader:
            try:
                x = float(row["POINT_X"])
                y = float(row["POINT_Y"])
            except:
                continue
            point = eps_scan.Point(x, y, False)
            points.append(point)
        max_cuttoff_net = 100000000000
        max_cuttoff_grid = 100000000000
        max_cuttoff_simple = 10000000000
        max_cuttoff_lin = 100000000000
        for j in xrange(args.repeats):
            rect = randomRect(points, args.r)

            dataset = []

            max_l_disc = 0
            total_red = 0
            total_blue = 0
            reg_red = 0
            reg_blue = 0
            for pt in points:
               if rect.inside(pt):
                   r_c = (args.red_bias + .5) > random.random()
                   reg_red += r_c
                   reg_blue += not r_c
               else:
                   r_c = .5 < random.random()
               total_red += r_c
               total_blue += not r_c
               dataset.append(pyscan.Point(pt.x, pt.y, r_c, not r_c))
            print total_red
            print total_blue
            """
            def inRegion(pt):
                return pt.getX() < -74.8  and 39.8 < pt.getY() and -76 < pt.getX() and pt.getY() < 40.15

            measured = [pt for pt in dataset if pt.getRedWeight() and inRegion(pt)]
            baseline = [pt for pt in dataset if not pt.getRedWeight() and inRegion(pt)]
            ms = measured #random.sample(measured, args.scatter)
            bs = baseline #random.sample(baseline, args.scatter)
            xa = [pt.getX() for pt in ms]
            xb = [pt.getX() for pt in bs]
            ya = [pt.getY() for pt in ms]
            yb = [pt.getY() for pt in bs]
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(xa, ya, marker=',', color='red')
            ax.scatter(xb, yb, marker=',', color='blue')
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            plt.show()
            """

            print reg_red + reg_blue
            print reg_red
            print reg_blue
            print total_blue - reg_blue
            print total_red - reg_red
            max_l_disc = (args.a * reg_red) / total_red + (args.b * reg_blue) / total_blue

            actualGridSub = []
            approxGridSub = []

            def maxD(rA, mx):
                return mx - rA.fValue()


            for i in xrange(1, args.repr + 1):
                k = int(2**i) * args.n
                print("k = " + str(k))

                grid = pyscan.makeGrid(dataset, k)
                if k < max_cuttoff_lin:
                    st = time.time()
                    print "Linear"
                    regionApprox = pyscan.maxSubgridLinear(grid, k * i, args.a, args.b)
                    et = time.time() - st
                    if et > args.max_time:
                        max_cuttoff_lin = k
                    else:
                        timingWriter.writerow({"time": et, "alg": "ff", "r": k,
                                               "maxDiff": maxSymDiff(rect, regionApprox, dataset, grid)})
                if i < max_cuttoff_simple:
                    st = time.time()
                    print "Linear Simple"
                    regionApprox = pyscan.maxSubgridLinearSimple(grid, args.a, args.b)
                    et = time.time() - st
                    if et > args.max_time:
                        max_cuttoff_simple = i
                    else:
                        timingWriter.writerow({"time": time.time() - st, "alg": "fs", "r": i * args.n,
                                               "maxDiff" : maxSymDiff(rect, regionApprox, dataset, grid)})

                if not args.disable_slow:
                    if i < 5:
                        st = time.time()
                        regionApprox = pyscan.maxSubgridLinearSlow(grid, args.a, args.b)
                        et = time.time() - st
                        if et > args.max_time:
                            max_cuttoff_grid = k
                        else:
                            timingWriter.writerow({"time": time.time() - st, "alg": "ss", "r": k,
                                                   "maxDiff" : maxSymDiff(rect, regionApprox, dataset, grid)})
                        print "Linear Slow"
                    if i < 5:
                        netGrid = pyscan.makeNetGrid(dataset, k)
                        print "Net Grid created"
                        st = time.time()
                        regionApprox = pyscan.maxSubgridLinearSlow(grid, args.a, args.b)
                        et = time.time() - st
                        if et > args.max_time:
                            max_cuttoff_net = k
                        else:
                            timingWriter.writerow({"time": time.time() - st, "alg": "ns",
                                               "r": k,
                                               "maxDiff": maxSymDiff(rect, regionApprox, dataset, grid)})
                        print "Linear Net"
        time_file.close()
