import random
import bisect
import math
import Partitioning as Part
import itertools
import csv
import time
import numpy as np
import test_set
import numpy.random as npr

from SeidelTree import Segment, Line, to_line
from collections import deque
from Cuttings import approx_eq


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
                u_b_l.append((u_s, bw, rw))
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


def order_function(p1, p2):
    x = p2[0] - p1[0]
    y = p2[1] - p1[1]
    if y >= 0:
        return math.atan2(y, x)
    else:
        return math.pi - math.atan2(y, x)


def line_error_test(n, red_points, red_weights, blue_points, blue_weights):
    """
    :param big_samp:
    :param small_samp:
    :param weights:
    :param pt_num:
    :param n: The sub-sampling size of the small sample.
    :return:
    """
    net_sample = random.sample(blue_points, n // 2) + random.sample(red_points, n // 2)

    for i in range(len(net_sample)):
        sample_part = net_sample[i + 1:]


        order_f = lambda x: order_function(net_sample[i], x)

        red_delta = [0] * (len(sample_part) + 1)
        blue_delta = [0] * (len(sample_part) + 1)

        sample_part.sort(key=lambda x: order_f(x))
        angles = [order_f(p) for p in sample_part]


        for p_1, w in zip(red_points, red_weights):
            insertion_pt = bisect.bisect_left(angles, order_f(p_1))
            red_delta[insertion_pt] += math.copysign(w, p_1[1])

        for p_1, w in zip(blue_points, blue_weights):
            insertion_pt = bisect.bisect_left(angles, order_f(p_1))
            blue_delta[insertion_pt] += math.copysign(w, p_1[1])

        big_sample_curr = sum(1 for p in red_points if p[1] <= 0)
        small_sample_curr = sum(1 for p in blue_points if p[1] <= 0)

        max_discrepancy = 0
        max_line = Line(0, 0)
        a = 1.0 / sum(red_weights)
        b = 1.0 / sum(blue_weights)
        for db, ds, p_2 in zip(red_delta, blue_delta, sample_part):
            big_sample_curr += db
            small_sample_curr += ds

            if max_discrepancy <= abs(big_sample_curr * a - small_sample_curr * b):
                max_line = to_line(net_sample[i], p_2)
                max_discrepancy = abs(big_sample_curr * a - small_sample_curr * b)

        return max_discrepancy, max_line


def naive_line_error_test(n, red_points, red_weights, blue_points, blue_weights):
    net_sample = random.sample(blue_points, n / 2) + random.sample(red_points, n / 2)
    def gen():
        for p1, p2 in itertools.combinations(net_sample, 2):
            l = to_line(p1, p2)
            r = sum(w for p, w in zip(red_points, red_weights) if l.pt_eq_below(p))
            b = sum(w for p, w in zip(blue_points, blue_weights) if l.pt_eq_below(p))
            yield abs(r - b), l
    return max(gen(), key= lambda x: x[0])



def line_discrepancy(n, red_points, red_weights, blue_points, blue_weights, scan_alg="fast"):
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


    output_file = "{}_{}_{}_{}_{}_{}_{}.csv".format(vparam, sampling_alg, scan_alg, c, input_size, l_s, h_s)


    fieldnames = ["vparam", "input_size", "sampling_alg", "scan_alg", "n", "s", "time", "m_disc_approx", "m_disc"]
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for i in np.linspace(l_s, h_s, count):
            if vparam == "eps":
                eps = i
            elif vparam == "c":
                c = i
            elif vparam == "input_size":
                input_size = i
            else:
                raise ValueError("vparam=%s"%(vparam,))

            red_input_set = random.sample(r_points, input_size)
            blue_input_set = random.sample(b_points, input_size)

            #(min_cell_size, len(input_set), r, cell_size)
            n = int(round(1 / eps) + .1)

            print(n, eps)
            start_time = time.time()

            s = int(round(c / eps ** 2) + .1)
            red_points_s = random.sample(red_input_set, min(s, len(red_input_set)))
            blue_points_s = random.sample(blue_input_set, min(s, len(blue_input_set)))

            if sampling_alg == "chan":
                s = int(round(c / (eps ** (4/3)) * math.log(1/ eps) ** (2/3)) + .1)
                print(len(red_points_s)/s)
                red_points_s, red_weights = Part.chan_partitions(red_points_s, 2, min_cell_size=len(red_points_s) / s,
                                                    test_set_f=test_set.test_set_points)
                blue_points_s, blue_weights = Part.chan_partitions(blue_points_s, 2, min_cell_size=len(blue_points_s) / s,
                                                     test_set_f=test_set.test_set_points)
            elif sampling_alg == "quad":
                expon = 1 / (2 - math.log(3, 4))
                s = int(round(c / eps ** (2 * expon) * math.log(1/eps) ** expon) + .1)
                red_points_s, red_weights = Part.quadTreeSample(red_points_s, len(red_points_s) / s)
                blue_points_s, blue_weights = Part.quadTreeSample(blue_points_s, len(blue_points_s) / s)
            elif sampling_alg == "naive":
                red_weights = [1.0] * len(red_points_s)
                blue_weights = [1.0] * len(blue_points_s)
            else:
                raise ValueError("sampling_alg = {}".format(sampling_alg))

            print(n, len(red_points_s), len(blue_points_s))
            mx, l = line_discrepancy(n, red_points_s, red_weights, blue_points_s, blue_weights, scan_alg=scan_alg)
            end_time = time.time()

            actual_mx = exact_discrepancy(l, red_points, red_weights, blue_points, blue_weights)

            row = {"vparam": vparam, "input_size":input_size, "sampling_alg":sampling_alg, "scan_alg":scan_alg,
                   "n":n, "s":s, "time":end_time - start_time, "m_disc_approx": mx, "m_disc": actual_mx}
            writer.writerow(row)
            print(row)


def exact_discrepancy(line, red_points, red_weights, blue_points, blue_weights):
    r = sum(w for p, w in zip(red_points, red_weights) if line.pt_eq_below(p))
    b = sum(w for p, w in zip(blue_points, blue_weights) if line.pt_eq_below(p))
    return abs(r / sum(red_weights) - b / sum(blue_weights))


def generate_blue_and_red(s_size, q, r, eps=.1):
    full_points = [(random.random(), random.random()) for _ in range(s_size)]

    while True:
        p1 = random.choice(full_points)
        p2 = random.choice(full_points)
        print("got here")
        try:
            l = to_line(p1, p2)
            below_line = sum(1 for p in full_points if l.pt_eq_below(p))
            print(abs(below_line / len(full_points)))
            if abs(below_line / len(full_points) - r) <= eps:
                break
        except ZeroDivisionError:
            continue

    red_points = []
    blue_points = []
    for p in full_points:
        if l.pt_eq_below(p):
            if .5 - q / 2 < random.random():
                red_points.append(p)
            else:
                blue_points.append(p)
        else:
            if .5 + q / 2 < random.random():
                red_points.append(p)
            else:
                blue_points.append(p)
    return red_points, blue_points


if __name__ == "__main__":
    print("running experiment")
    red_points, blue_points = generate_blue_and_red(10000000, .01, .5)
    # red_points = [(random.random(), random.random()) for _ in range(1000000)]
    # blue_points = [(random.random(), random.random()) for _ in range(1000000)]

    for s_alg in ["chan", "naive", "quad"]:
        for sc_alg in ["fast"]:
            print(s_alg + " " + sc_alg)
            testing_framework(red_points, blue_points, .001, .0001, 10, input_size=min(len(red_points), len(blue_points)), sampling_alg=s_alg, scan_alg=sc_alg)