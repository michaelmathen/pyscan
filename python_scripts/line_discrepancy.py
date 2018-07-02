import random
import bisect
import math
import partitioning as Part
from partitioning import line_discrepancy, stat, exact_discrepancy
import itertools
import csv
import time
import numpy as np
import test_set
import numpy.random as npr

from seidel_tree import Segment, Line, to_line
from collections import deque



def horizontal_split_vertices(points, segment):
    up = []
    down = []
    for p in points:
        if segment.pt_eq_below(p):
            down.append(p)
        else:
            up.append(p)
    return up, down


def score(end_point, w_lines):
    b_score = 0
    r_score = 0
    for line, rw, bw in w_lines:
        if line.pt_eq_below(end_point):
            r_score += rw
            b_score += bw
    return b_score, r_score


class Polygon:

    def __init__(self, w_lines, b_below= 0, r_below=0):
        self.w_lines = w_lines
        self.weight = sum(abs(rw) + abs(bw) for _, rw, bw in w_lines)
        self.blue_below = b_below
        self.red_below = r_below

    def horz_split(self, segment):

        lower_blue = 0
        lower_red = 0
        u_b_l = []
        l_b_l = []
        for l, bw, rw in self.w_lines:
            if segment.same_line(l):
                continue
            elif l.crossed_by(segment):
                u_s, l_s = l.simple_split(segment)
                if u_s is not None:
                    u_b_l.append((u_s, bw, rw))
                if l_s is not None:
                    l_b_l.append((l_s, bw, rw))
            elif l.above_closed(segment):
                u_b_l.append((l, bw, rw))
            else:
                lower_blue += bw
                lower_red += rw
                l_b_l.append((l, bw, rw))

        scored_pts = []
        if segment.xl != -math.inf:
            b_s, r_s = score(segment.left_vertex, self.w_lines)
            scored_pts.append((segment.left_vertex, b_s + self.blue_below, r_s + self.red_below))
        if segment.xr != math.inf:
            b_s, r_s = score(segment.right_vertex, self.w_lines)
            scored_pts.append((segment.right_vertex, b_s + self.blue_below, r_s + self.red_below))

        return Polygon(u_b_l, b_below=(self.blue_below + lower_blue), r_below=(self.red_below + lower_red)), \
                Polygon(l_b_l, b_below=self.blue_below, r_below=self.red_below), \
                scored_pts


    def find_pretty_good_split_l(self):
        vals = [abs(bw) + abs(rw) for _, bw, rw in self.w_lines]
        total_w = sum(vals)
        p = [w / total_w for w in vals]
        segments = npr.choice([l for l, _, _ in self.w_lines], p = p)
        return segments

    def get_weight(self) -> float:
        return self.weight


def discrepancy_fast(r, red_points, blue_points):
    """
    Computes the discrepancy by computing a cutting in the dual space of the points.
    :param self:
    :param r:
    :param red_points:
    :param blue_points:
    :return:
    """
    total_weight = len(red_points) + len(blue_points)
    min_weight = total_weight / r

    w_lines = []
    for pt in red_points:
        w_lines.append((Segment(Line(-pt[0], pt[1]), -math.inf, math.inf), 0, 1.0))
    for pt in blue_points:
        w_lines.append((Segment(Line(-pt[0], pt[1]), -math.inf, math.inf), 1.0, 0))

    node_stack = deque([Polygon(w_lines=w_lines)])

    max_disc = 0
    max_pt = (0, 0)
    while node_stack:
        curr_node = node_stack.pop()
        if curr_node.get_weight() > min_weight:
            segment = curr_node.find_pretty_good_split_l()
            upper, lower, new_lines = curr_node.horz_split(segment)

            for p in new_lines:
                if abs(p[1] - p[2]) > max_disc:
                    max_pt = p[0]
                    max_disc = abs(p[1] - p[2])

            node_stack.append(upper)
            node_stack.append(lower)

    return max_disc, Line(-max_pt[0], max_pt[1])




def bottom_quad(p1, p2):
    return p2[1] - p1[1]




def line_error_test(n, red_points, red_weights, blue_points, blue_weights):
    r_count = max(min(n // 2, len(red_points)), 1)
    b_count = max(min(n // 2, len(blue_points)), 1)
    r_t = sum(red_weights)
    b_t = sum(blue_weights)
    r_p = [w / r_t for w in red_weights]
    b_p = [w / b_t for w in blue_weights]
    net_1 = [blue_points[i] for i in np.random.choice(range(len(blue_weights)), b_count, p=b_p)]
    net_2 = [red_points[i] for i in np.random.choice(range(len(red_weights)), r_count, p=r_p)]

    def stat(m, b):
        return abs(m - b)

    return line_discrepancy(net_1 + net_2, red_points, red_weights, blue_points, blue_weights, disc=stat)


def naive_line_error_test(n, red_points, red_weights, blue_points, blue_weights):
    net_sample = random.sample(blue_points, n / 2) + random.sample(red_points, n / 2)
    nr = sum(red_weights)
    nb = sum(blue_weights)
    def gen():
        for p1, p2 in itertools.combinations(net_sample, 2):
            l = to_line(p1, p2)
            r = sum(w for p, w in zip(red_points, red_weights) if l.pt_eq_below(p))
            b = sum(w for p, w in zip(blue_points, blue_weights) if l.pt_eq_below(p))
            yield stat(r / nr, b / nb), l
    return max(gen(), key= lambda x: x[0])



def line_test(n, red_points, red_weights, blue_points, blue_weights, scan_alg="fast"):
    if scan_alg == "fast":
        return line_error_test(n, red_points, red_weights, blue_points, blue_weights)
    elif scan_alg == "naive":
        return naive_line_error_test(n, red_points, red_weights, blue_points, blue_weights)
    else:
        return ValueError("Bad scan algorithm")
    # else:
    #     return discrepancy_fast(n, red_points, blue_points)



def testing_framework(r_points, b_points, l_s, h_s, count,
                      vparam="eps",
                      eps=.01,
                      c=1,
                      input_size=100000,
                      sampling_alg="naive",
                      scan_alg="fast"):


    output_file = "{}_discrepancy.csv".format(sampling_alg)


    fieldnames = ["vparam", "input_size", "sampling_alg", "scan_alg", "n", "s", "time", "m_disc_approx", "m_disc"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for eps in np.logspace(l_s, h_s, count):

            #eps = i


            red_input_set = random.sample(r_points, input_size)
            blue_input_set = random.sample(b_points, input_size)

            #(min_cell_size, len(input_set), r, cell_size)
            #n = int(round(1 / eps) + .1)
            n = 400
            start_time = time.time()

            s = int(round(c / eps ** 2) + .1)
            red_points_s = random.sample(red_input_set, min(s, len(red_input_set)))
            blue_points_s = random.sample(blue_input_set, min(s, len(blue_input_set)))
            if sampling_alg == "chans":
                s = int(round(c / (eps ** (4/3)) * math.log(1/ eps) ** (2/3)) + .1)
                print(len(red_points_s)/s)

                red_tree = Part.chan_partitions_simple(red_points_s, 22, min_cell_size=len(red_points_s) / s)
                red_points_s, red_weights_s = Part.sample_cells([leaf.get_points() for leaf in red_tree.get_leaves()], 1)
                blue_tree = Part.chan_partitions_simple(blue_points_s, 22, min_cell_size=len(blue_points_s) / s)
                blue_points_s, blue_weights_s = Part.sample_cells([leaf.get_points() for leaf in blue_tree.get_leaves()], 1)

            elif sampling_alg == "chan":
                s = int(round(c / (eps ** (4 / 3)) * math.log(1 / eps) ** (2 / 3)) + .1)
                print(len(red_points_s) / s)

                red_tree = Part.chan_partitions2(red_points_s, 22, min_cell_size=len(red_points_s) / s)
                red_points_s, red_weights_s = Part.sample_cells([leaf.get_points() for leaf in red_tree.get_leaves()], 1)
                blue_tree = Part.chan_partitions2(blue_points_s, 22, min_cell_size=len(blue_points_s) / s)
                blue_points_s, blue_weights_s = Part.sample_cells([leaf.get_points() for leaf in blue_tree.get_leaves()], 1)
            elif sampling_alg == "mat":
                s = int(round(c / (eps ** (4 / 3)) * math.log(1 / eps) ** (2 / 3)) + .1)
                print(len(red_points_s) / s)

                red_points_s, red_weights_s = Part.partitions(red_points_s, 16,
                                                             min_cell_size=len(red_points_s) / s,
                                                            test_set_f=test_set.test_set_points)
                blue_points_s, blue_weights_s = Part.partitions(blue_points_s, 16,
                                                            min_cell_size=len(blue_points_s) / s,
                                                              test_set_f=test_set.test_set_points)
            elif sampling_alg == "quad":
                expon = 1 / (2 - .5)
                s = int(round(c / eps ** (2 * expon) * math.log(1/eps) ** expon) + .1)
                red_points_s, red_weights_s = Part.quadTreeSample(red_points_s, len(red_points_s) / s)
                blue_points_s, blue_weights_s = Part.quadTreeSample(blue_points_s, len(blue_points_s) / s)
            elif sampling_alg == "naive":
                red_weights_s = [1.0] * len(red_points_s)
                blue_weights_s = [1.0] * len(blue_points_s)
            else:
                raise ValueError("sampling_alg = {}".format(sampling_alg))

            r_count = max(min(n // 2, len(red_points)), 1)
            b_count = max(min(n // 2, len(blue_points)), 1)
            r_t = sum(red_weights_s)
            b_t = sum(blue_weights_s)
            r_p = [w / r_t for w in red_weights_s]
            b_p = [w / b_t for w in blue_weights_s]
            net_1 = [blue_points[i] for i in np.random.choice(range(len(blue_weights_s)), b_count, p=b_p)]
            net_2 = [red_points[i] for i in np.random.choice(range(len(red_weights_s)), r_count, p=r_p)]
            mx, mxline = line_discrepancy(net_1 + net_2, red_points_s, red_weights_s, blue_points_s, blue_weights_s, lambda m, b: abs(m - b))
            end_time = time.time()

            actual_mx = exact_discrepancy(mxline, red_input_set, [1 for p in red_input_set], blue_input_set, [1 for p in blue_input_set])

            row = {"vparam": vparam, "input_size":input_size, "sampling_alg":sampling_alg, "scan_alg":scan_alg,
                   "n":n, "s":s, "time":end_time - start_time, "m_disc_approx": mx, "m_disc": actual_mx}
            writer.writerow(row)
            print(row)





def line_planted_test(points, r):
    """
    :param big_samp:
    :param small_samp:
    :param weights:
    :param pt_num:
    :param n: The sub-sampling size of the small sample.
    :return:
    """
    net_sample = random.sample(points, int(2/r + 1))

    md, ml = line_discrepancy(net_sample, points, [1] * len(points), [], [], disc=lambda m, b: -abs(m - r))
    return ml


def generate_blue_and_red(full_points, q, r):

    # while True:
    #     p1 = random.choice(full_points)
    #     p2 = random.choice(full_points)
    #     print("got here")
    #     try:
    #         l = to_line(p1, p2)
    #         below_line = sum(1 for p in full_points if l.pt_eq_below(p))
    #         print(abs(below_line / len(full_points)))
    #         if abs(below_line / len(full_points) - r) <= eps:
    #             break
    #     except ZeroDivisionError:
    #         continue
    l = line_planted_test(full_points, r)

    red_points = []
    blue_points = []
    red_below = 0
    blue_below = 0
    below = 0
    for p in full_points:
        if l.pt_eq_below(p):
            below += 1
            if .5 + q < random.random():
                red_below += 1
                red_points.append(p)
            else:
                blue_below += 1
                blue_points.append(p)
        else:
            if .5 < random.random():
                red_points.append(p)
            else:
                blue_points.append(p)
    print(below / len(full_points))
    print(abs(red_below / len(red_points) - blue_below / len(blue_points)))
    return red_points, blue_points


if __name__ == "__main__":
    print("running experiment")
    import partitioning_tests
    pts = partitioning_tests.upload_crimes("crimes.csv")


    #pts = [(random.random(), random.random()) for _ in range(1000000)]
    red_points, blue_points = generate_blue_and_red(pts, .2, .02)

    # red_points = random.sample(red_points, 100000)
    # blue_points = random.sample(blue_points, 100000)

    # import matplotlib.pyplot as plt
    # f, ax = plt.subplots()
    # x, y = zip(*red_points)
    # ax.scatter(x, y, color="red")
    # x, y = zip(*blue_points)
    # ax.scatter(x, y, color="blue")
    #
    # plt.show()

    for s_alg in ["quad", "naive", "mat", "chans", "chan"]:
        for sc_alg in ["fast"]:
            print(s_alg + " " + sc_alg)
            testing_framework(red_points, blue_points, -.5, -3, 80, input_size=min(len(red_points), len(blue_points)), sampling_alg=s_alg, scan_alg=sc_alg)