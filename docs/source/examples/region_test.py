
import shapefile
import pyscan
import matplotlib.pyplot as plt
import csv



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
        
    
alpha = .02
r_min = .05


core_set_pts2010 = pyscan.polygon_sample(regions, weights2010, 1000)
core_set_pts2017 = pyscan.polygon_sample(regions, weights2017, 1000)

net = pyscan.my_sample(core_set_pts2017, 200) + pyscan.my_sample(core_set_pts2010, 200)
disc_f = pyscan.DISC
disk, d_val = pyscan.max_disk_scale(net, 
                                  [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2017], 
                                  [pyscan.WPoint(1.0, p[0], p[1], 1.0) for p in core_set_pts2010],
                                  1,
                                  disc_f)
_, ax = plt.subplots(figsize=(16, 12))
plt.axis('off')
plot_points(ax, core_set_pts2010, "r")
plot_points(ax, core_set_pts2017, "b")
d = plt.Circle(disk.get_origin(), disk.get_radius(), color='g', alpha=.8)
ax.add_artist(d)
print(disk)
plt.show()
