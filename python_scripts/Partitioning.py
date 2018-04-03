
from SeidelTree import Trapezoid, to_line, Segment, Line, \
    approx_above, approx_eq_above

import PolyTree as Pt
import random
from collections import deque, Counter
import test_set
import matplotlib.pyplot as plt
import math

import statistics
import itertools


from time import time

def sample_cells(cells, cell_sample_size):
    tw = sum(1 for _ in itertools.chain.from_iterable(cells))

    out_pts = []
    weights = []
    for cell in cells:
        if cell:
            c_size = min(len(cell), int(.5 + round(cell_sample_size)))
            out_pts.extend(random.sample(cell, c_size))
            weights.extend([len(cell) / (float(tw) * c_size)] * c_size)

    return list(out_pts), list(weights)

class PartTree:

    def __init__(self, trapezoid):
        self.trap = trapezoid
        self.children = []

    def insert_children(self, c):
        self.children.append(c)

    def visualize(self, ax, min_x, max_x, min_y, max_y, viz_lines=False, color='r', num_levels=None):
        if num_levels is None or num_levels > 0:
            if num_levels is not None:
                num_levels -= 1
            self.trap.visualize(ax, min_x, max_x, min_y, max_y, viz_lines=viz_lines, color=color)
            for c in self.children:
                c.visualize(ax, min_x, max_x, min_y, max_y, viz_lines=viz_lines, color=color, num_levels=num_levels)


def random_test_set(pts, t, c=10):
    test_set = []
    for i in range(int(t * t) * c):
        p1, p2 = random.sample(pts, 2)
        test_set.append(to_line(p1, p2))
    return test_set



def partitions(pts, t, min_cell_size=100, cell_sample_size=1, cutting_f=Pt.compute_cutting, test_set_f=random_test_set, c=6):
    #print(t)
    final_cells = []
    cell_queue = deque([pts])
    #t = max(int(c * math.sqrt(r)), 2)
    while cell_queue:

        curr_pts = cell_queue.pop()
        pts_not_in_cells = curr_pts[:]
        #print("divide cell %d, t=%d" % (len(curr_pts),t))
        test_set = test_set_f(pts_not_in_cells, t ** 2 * c)
        weight_map = {line: 1 for line in test_set}
        i = 0
        while len(pts_not_in_cells) > max(len(curr_pts) / (c * t ** 2), min_cell_size):
            i += 1
            # r_i = math.ceil(r / 2 ** i)
            t_i = max(t /  (math.sqrt(2)**i), 1)
            #print("r_i=%d, t_i=%d"%(r_i, t_i)t

            while len(pts_not_in_cells) > max(len(curr_pts) / (2**i), min_cell_size):
                #print(len(pts_not_in_cells), len(curr_pts), len(test_set))
                tree = cutting_f(test_set, weight_map, pts_not_in_cells, t_i)
                #print("Points left over %d, %d" % (len(pts_not_in_cells), len(curr_pts) / (2 ** i)))

                cell = tree.get_heaviest()
                #print("Number of points in %d"%(len(cell.get_points()),))
                for l in cell.get_lines():
                    weight_map[l] *= 2
                lkup = set(cell.get_points())
                pts_not_in_cells = [p for p in pts_not_in_cells if p not in lkup]
                if len(cell.get_points()) <= min_cell_size:
                    final_cells.append(cell.get_points())
                else:
                    cell_queue.append(cell.get_points())
        if len(pts_not_in_cells) <= min_cell_size:
            final_cells.append(pts_not_in_cells)
        else:
            cell_queue.append(pts_not_in_cells)

    return sample_cells(final_cells, cell_sample_size)



class BoxCell:

    def __init__(self, pts, lines, lx=-float("inf"), rx=float("inf"), by=-float("inf"), ty=float("inf")):
        self.lx = lx
        self.rx = rx
        self.ty = ty
        self.by = by
        self.pts = pts
        self.lines = lines

    def horz_cut(self, my):
        u_pts = [p for p in self.pts if my < p[1]]
        l_pts = [p for p in self.pts if p[1] <= my]
        u_lines = []
        l_lines = []
        sl = Segment(Line(0, my), self.lx, self.rx)
        for l in self.lines:
            if sl.crossed_by(l):
                u_lines.append(l)
                l_lines.append(l)
            elif sl.above_closed(l):
                l_lines.append(l)
            else:
                u_lines.append(l)
        return BoxCell(u_pts, u_lines, self.lx, self.rx, my, self.ty), \
             BoxCell(l_pts, l_lines, self.lx, self.rx, self.by, my)

    def vertical_cut(self, mx):
        r_pts = [p for p in self.pts if mx < p[0]]
        l_pts = [p for p in self.pts if p[0] <= mx]
        l_lines = []
        r_lines =[]
        for l in self.lines:
            y_val = l.evaluate(mx)
            if approx_above(self.by, y_val) and approx_above(y_val, self.ty):
                l_lines.append(l)
                r_lines.append(l)
            elif approx_eq_above(self.ty, y_val):
                if l.a <= 0:
                    r_lines.append(l)
                else:
                    l_lines.append(l)
            else:
                if l.a <= 0:
                    l_lines.append(l)
                else:
                    r_lines.append(l)
        return BoxCell(r_pts, r_lines, mx, self.rx, self.by, self.ty), \
            BoxCell(l_pts, l_lines, self.lx, mx, self.by, self.ty)

    def get_median(self, ix):
        ord_vals = map(lambda x: x[ix], self.pts)
        return statistics.median(ord_vals)

    def get_x_median(self):
        return self.get_median(0)

    def get_y_median(self):
        return self.get_median(1)

    def point_count(self):
        return len(self.pts)

    def get_points(self):
        return self.pts

    def get_lines(self):
        return self.lines


def point_cuts(pts, lines, max_number):
    """
    Divides the points into sets less than max_number.
    Uses alternating x, y order for the points. By doing this
    we get a set of cells with parameter n^7.9...
    :param pts:
    :param max_number:
    :return:
    """
    #Sort in the x-order
    out_cells = []
    cells = deque()
    cells.append((BoxCell(pts, lines), True))
    while cells:
        rect, order_x = cells.pop()
        if rect.point_count() <= max_number:
            out_cells.append(rect)
        else:
            if order_x:
                mv = rect.get_x_median()
                l_r, b_r = rect.vertical_cut(mv)
            else:
                mv = rect.get_y_median()
                l_r, b_r = rect.horz_cut(mv)
            cells.append((l_r, not order_x))
            cells.append((b_r, not order_x))
            #print(mv)
    return out_cells

"""
 TODO Move the poly tree algorithm over here and create a chan version of it.
 
 1) Need to modify the cutting so that the test sets are still internal to each cell. In other
 words the cutting needs to be constrained to the cell above it. (done)
 
 2) Change this so that we compute each cutting with a branching factor of b. This might involve
 doing a search over the cuttings.

3) Test set should be computed with size 1/eps^{4/3} log^{2/3 + 2} 1/eps. Or in other 
words the #cells * log^2(#cells). Or we could construct a sample of points of size 
sqrt(#cells) * log #cells and then feed this into the dual algorithm with r = sqrt(#cells) to get the lines.
This should scale like #cells * log #cells which is the optimal... 

"""

def chan_partitions(pts, r, min_cell_size=100, cell_sample_size=1,
                    cutting_f=Pt.compute_cutting,
                    test_set_f=test_set.test_set_dual_exact_t, c=1):

    final_cells = []

    n = (len(pts) / min_cell_size)
    test_set_size = n * (math.log(n) ** 2) * c
    print(test_set_size)
    test_set = test_set_f(pts, test_set_size)
    weight_map = {line: 1 for line in test_set}
    curr_level = [(pts, test_set)]

    t = 1
    b = 10
    while curr_level:
        t *= b

        active_level = []
        #inactive_cells = []

        for curr_pts, test_set in curr_level:
            #print(len(curr_pts) ,min_cell_size)
            if len(curr_pts) <= min_cell_size:
                final_cells.append(curr_pts)
            else:
                active_level.append((curr_pts, test_set))


        # ensure that the active level is processed in a random order
        random.shuffle(active_level)

        weight_map = {l: 1 for l in weight_map}
        curr_level = []
        for curr_pts, test_set in active_level:
            #print(len(curr_pts), len(test_set))

            sub_tree = cutting_f(test_set, weight_map, curr_pts, r)
            sub_partitions = sub_tree.get_leaves()

            #print("computed subtree %d"%(r,))
            sub_cells = []
            # Cut each node of the sub-tree into sub-cells containing at most 2n/(tb) points.
            #b = max(len(sub_partitions), 1)
            #b = r * r * 10
            part_size = max(2 * len(pts) / (t * b), min_cell_size)
            print("t = {} Psz = {} b = {} psz = {} lsz = {}".format(t, part_size, len(sub_partitions), len(curr_pts), len(test_set)))
            for trap in sub_partitions:
                cells_list = point_cuts(trap.get_points(), trap.get_lines(),  part_size)
                sub_cells.extend(cells_list)

            # Compute the count of line crossings
            l_counts = Counter()
            for sub_cell in sub_cells:
                l_counts.update(sub_cell.get_lines())

            # Apply re-weighting operation
            for l in l_counts:
                weight_map[l] *= (1 + 1.0 / len(sub_partitions)) ** (l_counts[l])

            for sub_cell in sub_cells:
                curr_level.append((sub_cell.get_points(), sub_cell.get_lines()))

    return sample_cells(final_cells, cell_sample_size)




def quadTreeSample(pts, min_cell_size=100, cell_sample_size=1):
    boxes = point_cuts(pts, [], min_cell_size)
    return sample_cells([b.get_points() for b in boxes], cell_sample_size)

#random_sample_time_plot2(999, 1000, 100000, 10)
#print("got here")
#time_plot2(200, 10000, 1000000, 20, r=2, c=.1)
#error_plot_output_size(100000, 100, 800, 10, r=2)
#random_sample_plot(100000, 200, 1200, 20)
#random_error_plot(10000, 20, 200, 20, r=2048, c=.1)
#random_error_plot2(10000, 200, 1, 10, 20, r=2048, c=.1)
#random_sample_time_plot(10000, 50, 1000, 100)
