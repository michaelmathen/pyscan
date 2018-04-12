
from SeidelTree import Trapezoid, to_line, Segment, Line, \
    approx_above, approx_eq_above

import PolyTree as Pt
import Cuttings as Ct
import random
from collections import deque, Counter

import math
import numpy.random as npr
import numpy as np
import test_set
import statistics
import itertools




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
        u_lines, l_lines= Pt.split_lines(sl, self.lines)
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
                e1 = Segment(l, l.xl, mx)
                e2 = Segment(l, mx, l.xr)
                l_lines.append(e1)
                r_lines.append(e2)
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


def chan_process_level(prev_level, part_size, test_set, cutting_f, r):
    # ensure that the active level is processed in a random order
    random.shuffle(prev_level)

    weight_map = {l: 1 for l in test_set}
    next_level = []
    for curr_pts, t_s in prev_level:
        # print(len(curr_pts), len(test_set))

        sub_tree = cutting_f(t_s, weight_map, curr_pts, r)
        sub_partitions = sub_tree.get_leaves()

        # print("computed subtree %d"%(r,))
        sub_cells = []
        # Cut each node of the sub-tree into sub-cells containing at most 2n/(tb) points.

        # print("t = {} Psz = {} b = {} psz = {} lsz = {}".format(t, part_size, len(sub_partitions), len(curr_pts), len(test_set)))
        for trap in sub_partitions:
            cells_list = point_cuts(trap.get_points(), trap.get_lines(), part_size)
            sub_cells.extend(cells_list)

        # Compute the count of line crossings
        l_counts = Counter()
        for sub_cell in sub_cells:
            l_counts.update(sub_cell.get_lines())

        # Apply re-weighting operation
        for l in l_counts:
            weight_map[l] *= (1 + 1.0 / len(sub_partitions)) ** (l_counts[l])

        for sub_cell in sub_cells:
            next_level.append((sub_cell.get_points(), sub_cell.get_lines()))

    return next_level


def find_light_cells(prev_level, min_cell_size):
    final_cells = []
    active_level = []
    for curr_pts, t_s in prev_level:
        if len(curr_pts) <= min_cell_size:
            final_cells.append(curr_pts)
        else:
            active_level.append((curr_pts, t_s))
    return  active_level, final_cells


# def chan_partitions(pts, r, min_cell_size=100, cell_sample_size=1,
#                     cutting_f=Pt.compute_cutting,
#                     test_set_f=test_set.test_set_dual_exact_t, c=1, c2=2, max_h=2):
#     b = r * r * c2
#     if max_h is None:
#         max_t = (len(pts) / min_cell_size)
#         max_h = len(pts)
#     else:
#         max_t = int(b ** (max_h) + .1)
#
#     active_levels = deque()
#     active_levels.append(pts)
#     final_cells = []
#
#     while active_levels:
#         curr_pts = active_levels.pop()
#         if len(curr_pts) <= min_cell_size:
#             final_cells.append(curr_pts)
#             continue
#         test_set_size = int(max_t * (math.log(max_t) ** 2) * c + .5)
#         print(test_set_size)
#         test_set = test_set_f(curr_pts, test_set_size)
#         curr_level = [(curr_pts, test_set)]
#         t = 1
#         for i in range(max_h):
#             max_part_size = max(2 * len(curr_pts) / (t * b), 1)
#             print(max_part_size)
#             curr_level, finished_cells = find_light_cells(curr_level, min_cell_size)
#             final_cells.extend(finished_cells)
#             if curr_level:
#                 curr_level = chan_process_level(curr_level, max_part_size, test_set, cutting_f, r)
#             else:
#                 break
#             t *= b
#         for cell_pts, _ in curr_level:
#             active_levels.append(cell_pts)
#     return sample_cells(final_cells, cell_sample_size)



def chan_partitions(pts, r, min_cell_size=100, cell_sample_size=1,
                    cutting_f=Pt.compute_cutting,
                    test_set_f=test_set.test_set_dual_exact_t, c=1):

    import test_set
    final_cells = []

    n = (len(pts) / min_cell_size)
    if test_set_f == test_set.test_set_dual_exact_t or test_set_f == test_set.test_set_dual:
        test_set_size = int(round(n) + .1)
    else:
        test_set_size = int(n * (math.log(n) ** 2) * c + .5)

    line_test_set = test_set_f(pts, test_set_size)
    print(len(line_test_set), int(n * (math.log(n) ** 2) * c + .5))
    weight_map = {}
    for l in line_test_set:
        weight_map[l] = 1
    curr_level = [(pts, line_test_set)]

    t = 1
    b = r * r
    #tc = 1
    while curr_level:


        active_level = []

        for curr_pts, t_s in curr_level:
            #print(len(curr_pts) ,min_cell_size)
            if len(curr_pts) <= min_cell_size:
                final_cells.append(curr_pts)
            else:
                active_level.append((curr_pts, t_s))


        # ensure that the active level is processed in a random order
        random.shuffle(active_level)
        t = len(final_cells) + len(active_level)

        for l in weight_map:
            weight_map[l] = 1

        curr_level = []
        #tc *= r
        #total_weight = sum(weight_map[l] for l in line_test_set)
        for curr_pts, t_s in active_level:
            #print(len(curr_pts), len(test_set))


            #sub_tree = Pt.compute_cutting_weight(t_s, weight_map, curr_pts, total_weight / tc)
            sub_tree = cutting_f(t_s, weight_map, curr_pts, r)
            sub_partitions = sub_tree.get_leaves()

            #print("computed subtree %d"%(r,))
            sub_cells = []
            # Cut each node of the sub-tree into sub-cells containing at most 2n/(tb) points.
            part_size = max(2 * len(pts) / (t * b), min_cell_size)
            #print("t = {} Psz = {} b = {} psz = {} lsz = {} wsz={}".format(t, part_size, len(sub_partitions), len(curr_pts), len(t_s), sum(weight_map[l] for l in t_s)))
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


def random_gap(p):
    u = max(random.random(0.0, 1.0), np.finfo(float).eps)
    return math.floor(math.log(u) / math.log(1-p))


def p_sample(data, p):
    i = random_gap(p, len(data))
    while i < len(data):
        yield data[i]
        i += random_gap(p)


def binom_rand_k(n, p):
    """
    This generates a random k from the binomial distribution.
    This takes time O(np) on expectation.
    :param n:
    :param p:
    :return:
    """
    log_q = math.log(1.0 - p)
    x = 0
    sum = 0
    while True:
        sum += math.log(random.random()) / (n - x)
        if sum < log_q:
            return x
        x += 1


def test_set_power(b, t, n, m):
    return -int(-round(math.log(math.sqrt(b) * math.log(n) / m, 2) + .5))


def cell_rate(b, t, n, m):
    return math.sqrt(b) * math.log(n) / math.sqrt(t)


def chan_partitions2(pts, b, min_cell_size=100, c=1):
    final_cells = []

    s = (len(pts) / min_cell_size)
    test_set_size = int(round(s) + .1)
    line_test_set, dual_tree = test_set.test_dual_tree(pts, test_set_size)
    m = len(line_test_set)
    n = len(pts)

    print(len(line_test_set), int(n * (math.log(n) ** 2) * c + .5))


    curr_test_set_power = test_set_power(b, 1, n, m)
    p = 2.0 ** curr_test_set_power
    r_hat = list(p_sample(line_test_set, p))
    lines_in_tree = set(r_hat)

    tree = Pt.PolyTree2(points=pts, lines=r_hat)
    curr_level = [(tree.root, None)]
    weight_map = {l:0 for l in line_test_set}
    r_weight_map = {l:0 for l in line_test_set}
    for l in r_hat:
        r_weight_map[l] = 1
    while curr_level:

        active_level = []
        for curr_root, p in curr_level:
            if len(curr_root.get_points()) <= min_cell_size:
                final_cells.append((curr_root, p))
            else:
                active_level.append((curr_root, p))

        # ensure that the active level is processed in a random order
        random.shuffle(active_level)

        #update the current working test set
        t = len(final_cells) + len(active_level)
        if test_set_power(b, t, n, m) != curr_test_set_power:
            curr_test_set_power = test_set_power(b, t, n, m)
            p = 2.0 ** curr_test_set_power

            for l in r_hat:
                r_weight_map[l] = 0
            new_r_hat = list(p_sample(line_test_set, 2.0 ** curr_test_set_power))
            for l in new_r_hat:
                r_weight_map[l] = 1

            lines_to_insert = set(new_r_hat) - lines_in_tree
            for l in lines_to_insert:
                tree.add_line(l)
            # add new lines
            lines_in_tree |= lines_to_insert

        curr_level = []
        for curr_root, parent in active_level:
            sub_partitions = tree.cutting_b(curr_root, parent, b / 4, weight_map)
            sub_cells = []
            # Cut each node of the sub-tree into sub-cells containing at most 2n/(tb) points.
            part_size = max(2 * len(pts) / (t * b), min_cell_size)
            #print("t = {} Psz = {} b = {} psz = {} lsz = {} wsz={}".format(t, part_size, len(sub_partitions), len(curr_pts), len(t_s), sum(weight_map[l] for l in t_s)))
            for cell, parent in sub_partitions:
                cell_list = tree.partition(cell, parent, part_size)
                sub_cells.extend(cell_list)

            # 3
            # Sample the cells of this decomposition and then count the lines inside of them.
            l_counts = Counter()
            for sub_cell, _ in p_sample(sub_cells, cell_rate(b, t, n, m)):
                # TODO implement these two methods
                dual_cell = sub_cell.to_dual()
                dual_pts = dual_tree.contained_vertices(dual_cell)

                new_lines = [Ct.to_dual_line(p) for p in dual_pts]
                l_counts.update(new_lines)

            # 3(a) and 3(b)
            for l in l_counts:
                old_weight = weight_map[l]
                weight_map[l] *= (1 + 1.0 / b) ** (l_counts[l])

                #simulates taking a sample from the newly added lines
                k = binom_rand_k(weight_map[l] - old_weight, p)
                if k > 0:
                    # add these newly sampled lines to the r_weight_map
                    r_weight_map[l] += k
                    if l not in lines_in_tree:
                        tree.add_line(l)
            curr_level.extend(sub_cells)

    return tree



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
