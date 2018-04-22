import pyscan
import argparse
import random
import numpy.random as npr
import csv
import bisect
import time
import eps_scan
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter


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
                        default=20,
                        help='The max grid approximation.')
    parser.add_argument('--save',
                        help='Where to save this.')
    parser.add_argument('--disable_slow',
                        action='store_true',
                        help="Disable the slow algorithm")
    parser.add_argument('--timing', action='store_true', help='Timing step')
    parser.add_argument('--scatter',
                        type=int,
                        default=1000,
                        help='The size of the random sample to use with the visualization')

    parser.add_argument('--prefix',
                        required=True,
                        help="prefix for the timing and error outputs")
    
    parser.add_argument('--query',
                        required=True,
                        help='Any crime with this query text will be labeled.')

    args = parser.parse_args()
    csv_points = []
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
            el_id = args.query in row["TEXT_GENERAL_CODE"].lower()
            point = pyscan.Point(x, y, el_id, not el_id)
            numAnomalies += el_id
            points.append(point)
            csv_points.append((x, y, row["TEXT_GENERAL_CODE"]))

        print("Number of measured: " + str(numAnomalies))
        print("Number of Points: " + str(len(points)))
        print("Percent: " + str(float(numAnomalies) / len(points)))


        actualGridSub = []
        approxGridSub = []
        if args.timing:
            fieldnames = ["time", "alg", "r"]
            time_file = open(args.prefix + "_time.csv", 'w')
            timingWriter = csv.DictWriter(time_file, fieldnames=fieldnames)
            timingWriter.writeheader()
            for i in xrange(1, args.repr + 1):
                print i
                grid = pyscan.makeGrid(points, args.n * i)
                st = time.time()
                regionApprox = pyscan.maxSubgridLinear(grid, i * args.n, args.a, args.b)
                timingWriter.writerow({"time": time.time() - st, "alg": "ff", "r": i * args.n})
                st = time.time()
                regionApprox = pyscan.maxSubgridLinearSimple(grid, args.a, args.b)
                timingWriter.writerow({"time": time.time() - st, "alg": "fs", "r": i * args.n})
                if not args.disable_slow:
                    st = time.time()
                    regionApprox = pyscan.maxSubgridLinearSlow(grid, args.a, args.b)
                    timingWriter.writerow({"time": time.time() - st, "alg": "ss", "r": i * args.n})

            time_file.close()
        else:
            fieldnames = ["error", 'r_prime', 'r']
            error_file = open(args.prefix + "_error.csv", 'w')
            errorWriter = csv.DictWriter(error_file, fieldnames=fieldnames)
            errorWriter.writeheader()
            points = random.sample(points, 2500)
            for i in xrange(1, args.repr + 1):
                print i
                error_sum = 0
                for j in xrange(50):
                    dataset = []
                    for pt in points:
                        r_c = .5 > random.random()
                        dataset.append(pyscan.Point(pt.getX(), pt.getY(), r_c, not r_c))
                    (a, b) = randomLinFunc()
                    grid = pyscan.makeGrid(dataset, args.n)

                    regionApprox1 = pyscan.maxSubgridLinear(grid, i * args.n / 2, a, b)
                    regionApprox2 = pyscan.maxSubgridLinearSimple(grid, a, b)
                    error_sum += abs(regionApprox1.fValue() - regionApprox2.fValue())
                errorWriter.writerow({"error": error_sum / 50.0,
                                      "r_prime": i * args.n / 2,
                                      "r": i})
            error_file.close()


        """
        for pt in csv_points:
            if region.inside(eps_scan.Point(pt[0], pt[1], True)):
                in_region_cases[pt[2]] += 1
            else:
                out_region_cases[pt[2]] += 1

        most_ratio = Counter()
        least_ratio = Counter()        
        for case in in_region_cases:
            most_ratio[case] = float(in_region_cases[case]) / (out_region_cases[case] + in_region_cases[case])
            least_ratio[case] = float(out_region_cases[case]) / (out_region_cases[case] + in_region_cases[case])            
        print(most_ratio.most_common())
        print(least_ratio.most_common())
        """
        """
        print region
        def inRegion(pt):
            #print(pt.getX() < -74.8  and 39.8 < pt.getY() and -76 < pt.getX() and pt.getY() < 40.15)
            return pt.getX() < -74.8  and 39.8 < pt.getY() and -76 < pt.getX() and pt.getY() < 40.15
        #pt.x

        measured = [pt for pt in points if pt.getRedWeight() and inRegion(pt)]
        baseline = [pt for pt in points if not pt.getRedWeight() and inRegion(pt)]
        ms = random.sample(measured, args.scatter)
        bs = random.sample(baseline, args.scatter)
        xa = [pt.getX() for pt in ms]
        xb = [pt.getX() for pt in bs]
        ya = [pt.getY() for pt in ms]
        yb = [pt.getY() for pt in bs]
        
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.scatter(xa, ya, marker='.', color='red')
        ax.scatter(xb, yb, marker='.', color='blue')
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])        
        #plt.gca().axison = False
        ax.add_patch(
                patches.Rectangle(
                    (region.lowCol(), region.lowRow()),   # (x,y)
                    region.upCol() - region.lowCol(),          # width
                    region.upRow() - region.lowRow(),
                    fill=False,
                    alpha=.5,
                    color="black",
                    linewidth=3)          # height
            )

        ax.axis("tight")
        if not (args.save is None):
            pyplot.savefig(args.save, bbox_inches='tight', frameon=False)
        else:
            pyplot.show()
        """