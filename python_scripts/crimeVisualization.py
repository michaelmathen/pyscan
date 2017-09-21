import pyscan
import argparse
import random
import csv
import bisect
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
from collections import Counter


"""
Given a crime keyword find a region that is interesting and display that region. 
Also explain how interesting that region is.
"""

def permSampleL(collection, s, n):
    s_a = random.sample(collection, s)
    s_b = random.sample(collection, s)
    N = random.sample(collection, n)
    return (N, s_a, s_b)

def mkCDF(statistic):
    values = sorted(statistic)
    return lambda x: float(bisect.bisect_left(values, x) - 1) / len(values)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--n', type=int,
                        default=100,
                        help='The size of the random net')
    parser.add_argument('--repr', type=int,
                        default=20,
                        help='The number of replications to use for the power tests.')

    parser.add_argument('--alg', 
                        default="linRect",
                        choices=["linRect", "linRectS"],
                        help='The name of the algorithm to use')
    parser.add_argument('--save',
                        help='Where to save this.')

    parser.add_argument('--scatter',
                        type=int,
                        default=1000,
                        help='The size of the random sample to use with the visualization')

    
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
        grid = pyscan.makeGrid(points, args.n)
        if args.alg == "linRect":
            region = pyscan.maxSubgridLinear(grid, args.n, -.2, .1)
        else:
            region = pyscan.maxSubgridLinearSimple(grid, -.2, .1)
        region = grid.toRectangle(region)
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
