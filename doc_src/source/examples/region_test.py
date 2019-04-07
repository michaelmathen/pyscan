import shapefile
import pyscan
import matplotlib.pyplot as plt
import csv
import itertools


def plot_points(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0] )
        ys.append(pt[1])
    ax.scatter(xs, ys, color=c, marker='.')

def plot_points_traj(ax, pts, c):
    xs = []
    ys = []
    for pt in pts:
        xs.append(pt[0])
        ys.append(pt[1])
    ax.plot(xs, ys, color=c)
    
def plot_approx(ax, regions, core_set_pts):
    for reg in regions:
        plot_points_traj(ax, reg, "g")
    plot_points(ax, core_set_pts, "b")
    ax.set_axis_off()


shape = shapefile.Reader("county_shapes/cb_2017_us_county_500k.shp")
population2017 = {}
population2010 = {}

with open("county_population/PEP_2017_PEPANNRES_with_ann.csv", encoding='latin-1') as f:
    reader = csv.DictReader(f)
    headers = next(reader, None)
    for row in reader:
        population2017[row['GEO.id2'][-3:]] = int(row['respop72017'])
        population2010[row['GEO.id2'][-3:]] = int(row['respop72010'])
        
regions = []
weights2017 = []
weights2010 = []
for reg in shape.shapeRecords():
    ignore = False
    for p in reg.shape.points:
        # remove counties outside of the continental US
        if not (-124.84 <= p[0] <= -66.9 and 24.396 <= p[1] <= 49.4):
            ignore = True
            break
    if not ignore:
        weights2010.append(population2010[reg.record[1]]) #reg.record[2], reg.record[5])
        weights2017.append(population2017[reg.record[1]]) #reg.record[2], reg.record[5])
        regions.append([pyscan.Point(p[0], p[1], 1.0) for p in reg.shape.points])
        
disc_f = pyscan.DISC
alpha = 0.5
r_min = 1.0
r_max = 4.0


disk, value = pyscan.max_disk_region(regions, regions, weights2010, regions, weights2017, r_min, r_max, alpha, disc_f)
print(disk)
