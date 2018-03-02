
from SeidelTree import compute_cutting, dual_cutting, Trapezoid, to_line
import random
from collections import deque
import matplotlib.pyplot as plt
import math
import numpy as np
from time import time

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



def partitions(pts, r, c, min_cell_size=100, cell_sample_size=1, test_set_f=random_test_set):
    #print(t)
    final_cells = []
    cell_queue = deque([pts])
    t = max(int(c * math.sqrt(r)), 2)
    while cell_queue:

        curr_pts = cell_queue.pop()
        pts_not_in_cells = curr_pts[:]
        #print("divide cell %d, t=%d" % (len(curr_pts),t))
        test_set = test_set_f(pts_not_in_cells, t)
        weight_map = {line: 1 for line in test_set}
        i = 0
        while len(pts_not_in_cells) > max(len(curr_pts) / r, min_cell_size):
            i += 1
            r_i = math.ceil(r / 2 ** i)
            t_i = max(int(c * math.sqrt(r_i)), 2)
            #print("r_i=%d, t_i=%d"%(r_i, t_i))

            while len(pts_not_in_cells) > max(len(curr_pts) / (2**i), min_cell_size):

                random.shuffle(test_set)
                tree = compute_cutting(test_set, weight_map, pts_not_in_cells, t_i)
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

    out_pts = []
    weights = []
    for cell in final_cells:
        if cell:
            c_size = min(len(cell), int(.5 + round(cell_sample_size)))
            out_pts.extend(random.sample(cell, c_size))
            weights.extend([len(cell) / (float(len(pts)) * c_size)] * c_size)

    return list(out_pts), list(weights)


def chan_partitions(pts, r, c, min_cell_size=100, cell_sample_size=1, test_set_f=random_test_set):
    """
    Computes a partitioning based on the scheme from CHAN10
    :param pts:
    :param r:
    :param c:
    :param min_cell_size:
    :param cell_sample_size:
    :param test_set_f:
    :return:
    """
    #print(t)
    final_cells = []
    cell_queue = deque([pts])
    t = max(int(c * math.sqrt(r)), 2)


    while cell_queue:

        curr_pts = cell_queue.pop()
        pts_not_in_cells = curr_pts[:]
        #print("divide cell %d, t=%d" % (len(curr_pts),t))
        test_set = test_set_f(pts_not_in_cells, t)
        weight_map = {line: 1 for line in test_set}
        tree = compute_cutting(test_set, weight_map, curr_pts, t)
        for t in tree.get_trapezoids():
            sub_tree = compute_cutting(t.get_lines(), weight_map, t.get_points(), b)

        i = 0
        while len(pts_not_in_cells) > max(len(curr_pts) / r, min_cell_size):
            i += 1
            r_i = math.ceil(r / 2 ** i)
            t_i = max(int(c * math.sqrt(r_i)), 2)
            #print("r_i=%d, t_i=%d"%(r_i, t_i))

            while len(pts_not_in_cells) > max(len(curr_pts) / (2**i), min_cell_size):

                random.shuffle(test_set)
                tree = compute_cutting(test_set, weight_map, pts_not_in_cells, t_i)
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

    out_pts = []
    weights = []
    for cell in final_cells:
        if cell:
            c_size = min(len(cell), int(.5 + round(cell_sample_size)))
            out_pts.extend(random.sample(cell, c_size))
            weights.extend([len(cell) / (float(len(pts)) * c_size)] * c_size)

    return list(out_pts), list(weights)


def line_error_test(big_samp, small_samp, weights, pt_num):
    test_pool = big_samp + small_samp
    for i in range(pt_num):
        [p1, p2] = random.sample(test_pool, 2)
        try:
            line = to_line(p1, p2)
            r_1 = sum(line.pt_eq_below(p) for p in big_samp) / float(len(big_samp))
            r_2 = sum(weights[i] for i, p in enumerate(small_samp) if line.pt_eq_below(p)) / sum(weights[i] for i, p in enumerate(small_samp))
        except ZeroDivisionError:
            continue

        yield abs(r_1 - r_2)


def line_error_test_unweighted(big_samp, small_samp, pt_num):
    test_pool = big_samp + small_samp
    for i in range(pt_num):
        [p1, p2] = random.sample(test_pool, 2)
        try:
            line = to_line(p1, p2)
            r_1 = sum(line.pt_eq_below(p) for p in big_samp) / float(len(big_samp))
            r_2 = sum(line.pt_eq_below(p) for p in small_samp) / float(len(small_samp))
        except ZeroDivisionError:
            continue

        yield abs(r_1 - r_2)


def random_sample_plot(point_count, l_s, h_s, k):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    errors = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        s_size = int(i + .5)
        sampled_pts = random.sample(big_sample, s_size)
        #print(sampled_weights)
        max_error = 0
        for e in line_error_test_unweighted(big_sample, sampled_pts,  100):
            max_error = max(max_error, e)
        print("finished run")
        errors.append(max_error)
        sizes.append(len(sampled_pts))
    f, a = plt.subplots()
    a.scatter(sizes, errors)
    plt.show()

def random_sample_time_plot(point_count, l_s, h_s, k):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    times = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        s_size = int(i + .5)
        s = time()
        sampled_pts = random.sample(big_sample, s_size)
        e = time()
        times.append(e - s)
        sizes.append(len(sampled_pts))
    f, a = plt.subplots()
    a.scatter(sizes, times)
    plt.show()

def random_sample_time_plot2(point_count, l_s, h_s, k):
    times = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(int(i))]
        s = time()
        sampled_pts = random.sample(big_sample, point_count)
        e = time()
        times.append(e - s)
        sizes.append(i)
    f, a = plt.subplots()
    a.scatter(sizes, times)
    plt.show()


def error_plot(point_count, l_s, h_s, k, r=4, c=2):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    errors = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=i)
        #print(sampled_weights)
        max_error = 0
        for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
            max_error = max(max_error, e)
        print("finished run %f"%(i, ))
        errors.append(max_error)
        sizes.append(len(sampled_pts))
    #print(sizes, errors)
    f, a = plt.subplots()
    a.scatter(sizes, errors)
    plt.show()


def time_plot(point_count, min_cell_size, l_s, h_s, k, r=4, c=2):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    times = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        s = time()
        sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=min_cell_size, cell_sample_size=i,
                                                  test_set_f=random_test_set)
        e = time()
        #print(sampled_weights)
        print("finished run %f"%(i, ))
        times.append(e - s)
        sizes.append(len(sampled_pts))
    #print(sizes, errors)
    f, a = plt.subplots()
    a.scatter(sizes, times)
    plt.show()

def time_plot2(output_size, l_s, h_s, k, r=4, c=2):
    times = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(int(i))]
        min_cell_size = i / output_size
        s = time()
        sampled_pts, sampled_weights = partitions(big_sample, r, c,
                                                  min_cell_size=min_cell_size,
                                                  cell_sample_size=1,
                                                  test_set_f=random_test_set)
        e = time()
        #print(sampled_weights)
        print("finished run %f"%(i, ))
        times.append(e - s)
        sizes.append(i)
    #print(sizes, errors)
    f, a = plt.subplots()
    a.scatter(sizes, times)
    plt.show()

def error_plot2(point_count, min_cell_size, l_s, h_s, k, r=4, c=2):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    errors = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=min_cell_size, cell_sample_size=i)
        #print(sampled_weights)
        max_error = 0
        for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
            max_error = max(max_error, e)
        print("finished run %f"%(i, ))
        errors.append(max_error)
        sizes.append(len(sampled_pts))
    #print(sizes, errors)
    f, a = plt.subplots()
    a.scatter(sizes, errors)
    plt.show()




def random_error_plot(point_count, l_s, h_s, k, r=4, c=2):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    errors = []
    sizes = []

    for i in np.linspace(l_s, h_s, k):
        sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=i, test_set_f=random_test_set)
        #print(sampled_weights)
        max_error = 0
        for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
            max_error = max(max_error, e)
        print("finished run %f"%(i, ))
        errors.append(max_error)
        sizes.append(len(sampled_pts))
    #print(sizes, errors)
    f, a = plt.subplots()
    a.scatter(sizes, errors)
    plt.show()


def random_error_plot2(point_count, min_cell_size, l_s, h_s, k, r=4, c=2):
    big_sample = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    errors = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        sampled_pts, sampled_weights = partitions(big_sample, r, c, min_cell_size=min_cell_size, cell_sample_size=i,
                                                  test_set_f=random_test_set)
        #print(sampled_weights)
        max_error = 0
        for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
            max_error = max(max_error, e)
        print("finished run %f"%(i, ))
        errors.append(max_error)
        sizes.append(len(sampled_pts))
    #print(sizes, errors)
    f, a = plt.subplots()
    a.scatter(sizes, errors)
    plt.show()

#random_sample_time_plot2(999, 1000, 100000, 10)
time_plot2(200, 10000, 100000, 8, r=2048, c=.1)
#random_error_plot(10000, 20, 200, 20, r=2048, c=.1)
#random_error_plot2(10000, 200, 1, 10, 20, r=2048, c=.1)
#random_sample_time_plot(10000, 50, 1000, 100)
