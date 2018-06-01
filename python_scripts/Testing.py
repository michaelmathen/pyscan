import pyscan
import random
import time

def toPoints(point_loc, red_blue):
    points = []
    for (x,y), rb in zip(point_loc, red_blue):
        r = rb
        b = 1 - rb
        points.append(pyscan.Point(x, y, r, b))
    return points

def actualSubgridValue(grid, subgrid, a, b):
    gridVal = 0

    for i in xrange(subgrid.upCol(), subgrid.lowCol() + 1):
        for j in xrange(subgrid.upRow(), subgrid.lowRow() + 1):
            gridVal += a * grid.redWeight(i, j) + b * grid.blueWeight(i, j)
    return gridVal


if __name__ == "__main__":

    point_count = 100000
    point_loc = [(random.random(), random.random()) for i in range(point_count)]
    red_blue = [random.randint(0, 1) for i in range(point_count)]
    points = toPoints(point_loc, red_blue)
    #print(points)
    n = 200
    grid = pyscan.makeGrid(points, n)
    grid = pyscan.makeNetGrid(points, n)
    #tree = pyscan.SlabTree(grid, n * n * n)

#    for i in xrange(n):
#        for j in xrange(n):
#     #       print str(-.1 * grid.redWeight(i, j) + .1 * grid.blueWeight(i, j)) + " ",
#        print

    #print tree.measure(0, 0, 2, 2, 1, 1);

    #print(pyscan.maxSubgridKullSlow(grid, .01))
    #st = time.time()

    #print(time.time() - st)
    #st = time.time()
    sl = pyscan.maxSubgridLinear(grid, n * n * n, -.1, .1)
    ss = pyscan.maxSubgridLinearSimple(grid, -.1, .1)

    print(str(actualSubgridValue(grid, sl, -.1, .1)) + " " + str(sl))
    print(str(actualSubgridValue(grid, ss, -.1, .1)) + " " + str(ss))


    #print(time.time() - st)
    #bloom = pyscan.BloomFilter(10, .05)
    #for i in range(100):
    #    for j in range(100):
    #        print(grid.redWeight(i, j))
    #        print(grid.blueWeight(i, j))

