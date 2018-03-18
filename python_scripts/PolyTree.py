from Cuttings import approx_above, \
    approx_eq, \
    approx_eq_above, \
    weighted_shuffle, \
    Line
import pprint
import itertools
import random
from collections import deque
from SeidelTree import Segment, to_line
import math
import pydot
import numpy as np
from typing import List, Tuple
import matplotlib.pyplot as plt

id = 0


class Node:
    def is_terminal(self):
        return False


def det2(a1, a2, b1, b2):
    return a1 * b2 - a2 * b1

#
# class HLine:
#
#     def __init__(self, a, b, c = 1):
#         self.a = a
#         self.b = b
#         self.c = c
#
#     def join(self, other):
#         a = det2(self.b, self.c, other.b, other.c)
#         b = det2(self.a, self.c, other.a, other.c)
#         c = det2(self.a, self.b, other.a, other.b)
#         return HLine(a, b, c)
#
#     def evaluate(self, other):
#         return HLine(self.a * other.a, self.b * other.b, self.c * other.c)
#
#     def x_intercept(self, other):
#         pt = self.join(other)
#         return pt.a / pt.c
#
#     # These compare the first or second coordinates of the Hline
#     def c1_lt(self, other) -> bool:
#         return self.a * other.c < other.a * self.c
#
#     def c2_lt(self, other) -> bool:
#         return self.b * other.c < other.b * self.c
#
#     def c1_lte(self, other) -> bool:
#         return self.a * other.c <= other.a * self.c
#
#     def c2_lte(self, other) -> bool:
#         return self.b * other.c <= other.b * self.c
#
#     def c1_approx_lt(self, other) -> bool:
#         return approx_above(self.a * other.c < other.a * self.c)
#
#     def c2_approx_lt(self, other) -> bool:
#         return approx_above(self.b * other.c, other.b * self.c)
#
#     def c1_approx_lte(self, other) -> bool:
#         return approx_eq_above(self.a * other.c < other.a * self.c)
#
#     def c2_approx_lte(self, other) -> bool:
#         return approx_eq_above(self.b * other.c, other.b * self.c)
#
#     # These test to see whether the dual is above or below the line.
#     def dual_lt(self, dual) -> bool:
#         return 0 < dual.a * self.a + dual.b * self.b + dual.c * self.c
#
#     def dual_lte(self, dual) -> bool:
#         return 0 <= dual.a * self.a + dual.b * self.b + dual.c * self.c
#
#     def dual_approx_lt(self, dual) -> bool:
#         return approx_above(0, dual.a * self.a + dual.b * self.b + dual.c * self.c)
#
#     def dual_approx_lte(self, dual) -> bool:
#         return approx_eq_above(0, dual.a * self.a + dual.b * self.b + dual.c * self.c)
#
#     # These test to see if a line is below a segment or crosses it.
#
#     def crossing_interval(self, other, xl, cl, xr, cr):
#         h_object = self.join(other)
#         return approx_above(xl * h_object.c, h_object.a * cl) and \
#                approx_above(h_object.a * cr, xr * h_object.c)
#
#     def crossing_interval_closed(self, other, xl, cl, xr, cr):
#         h_object = self.join(other)
#         return approx_eq_above(xl * h_object.c, h_object.a * cl) and \
#                approx_eq_above(h_object.a * cr, xr * h_object.c)
#
#
#     def above_closed_interval(self, other, xl, xr, cl=1, cr=1):
#         h_object = self.join(other)
#         if approx_eq(h_object.c, 0.0):
#             return approx_eq_above(other.c / other.b, self.c / self.b)
#         else:
#             br = approx_above(h_object.a / h_object.c, xr / cr)
#             bl = approx_above(xl / cl, h_object.a / h_object.c)
#             if bl and br:
#                 return False
#             elif approx_eq_above(xr / cr, h_object.a / h_object.c):
#                 return h_object.c < 0
#             else:
#                 return h_object.c > 0
#
#     def above_interval(self, other, xl, xr, cl=1, cr=1):
#         h_object = self.join(other)
#         if approx_eq(h_object.c, 0.0):
#             return approx_above(other.c / other.b, self.c / self.b)
#         else:
#             vr = det2(h_object.a, h_object.c, xr, cr)
#             vl = det2(h_object.a, h_object.c, xl, cl)
#
#             if approx_eq_above(vr, 0):
#                 return False
#             elif approx_eq_above(xr / cr, h_object.a / h_object.c):
#                 return h_object.c < 0
#             else:
#                 return h_object.c > 0
#
# class HSegment(HLine):
#     def __init__(self, a, b, c, xl, xr):
#         super(HSegment, self).__init__(a, b, c)
#         self.xl = xl
#         self.xr = xr

#


class Segment_Node(Node):

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, segment, up=None, down=None):
        self.up = up
        self.down = down
        self.segment = segment

    def segment_b(self, seg):
        return seg.below_closed_interval(self.segment, seg.xl, seg.xr)

    def segment_a(self, seg):
        return seg.above_closed_interval(self.segment, seg.xl, seg.xr)

    def get_a(self):
        return self.up

    def get_b(self):
        return self.down

    def set_a(self, u):
        self.up = u

    def set_b(self, b):
        self.down = b

    def crosses(self, seg):
        """
        Has to intersect the interior of the segment.
        :param seg:
        :return:
        """
        if approx_eq(seg.a, self.segment.a):
            return False
        x_v = seg.x_intercept(self.segment)
        return approx_above(self.segment.xl, x_v) and \
               approx_above(x_v, self.segment.xr) and \
               approx_above(seg.xl, x_v) and \
               approx_above(x_v, seg.xr)



def horzontal_split_lines(lines, segment, key=lambda x: x):

    u_l = []
    l_l = []
    for l in lines:
        if segment.crossed_by(key(l)):
            u_l.append(l)
            l_l.append(l)
        else:
            if segment.same_line(key(l)):
                continue
            elif segment.above_closed(key(l)):
                l_l.append(l)
            elif segment.below_closed(key(l)):
                u_l.append(l)
    return u_l, l_l

def horizontal_split_vertices(points, segment):
    pt_below = []
    pt_above = []

    for p in points:
        if segment.pt_eq_below(p):
            pt_below.append(p)
        if segment.pt_eq_above(p):
            pt_above.append(p)

    return pt_above, pt_below

poly_id = 0
class Polygon(Node):
    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, border_lines=list(), w_lines=list(), points=list(), k=1000):
        self.border_lines = border_lines
        self.w_lines = w_lines
        self.points = points
        self.weight = sum(w for _, w in w_lines)
        global poly_id
        self.id = poly_id
        poly_id += 1
        self.k = k

    def visualize(self, ax, min_x, max_x, min_y, max_y):
        xs, ys = zip(*self.points)
        ax.scatter(xs, ys)
        tmp_border = [self.border_lines[-1]] + self.border_lines + [self.border_lines[0]]

        vertices_in_order = [ll.crossing_pt(lx) for ll, lx in zip(tmp_border[:-1], tmp_border[1:])]
        #zip(vertices_in_order)

    def get_border_vertices(self):
        tmp_border = [self.border_lines[-1]] + self.border_lines + [self.border_lines[0]]
        return [ll.crossing_pt(lx) for ll, lx in zip(tmp_border[:-1], tmp_border[1:])]

    def get_points(self):
        return self.points

    def get_lines(self):
        return [l for l, _ in self.w_lines]

    def pt_count(self):
        return len(self.points)

    def is_terminal(self):
        return True

    def horz_split(self, segment):
        """
        This can only cut in two, because we always do vertical splits first.
        :param segment:
        :return:
        """
        up = []
        down = []
        for p in self.points:
            if segment.pt_eq_below(p):
                down.append(p)
            else:
                up.append(p)
        u_l = []
        l_l = []
        #print(self.w_lines)
        for l, w in self.w_lines:
            if segment.crossed_by(l):
                u_s, l_s = l.simple_split(segment)
                if u_s is not None: u_l.append((u_s, w))
                if l_l is not None: l_l.append((l_s, w))
            else:
                if segment.same_line(l):
                    continue
                elif l.above_closed(segment):
                    u_l.append((l, w))
                else:
                    l_l.append((l, w))

        # preserves the order of the lines
        upper_border_lines, lower_border_lines = horzontal_split_lines(self.border_lines, segment)
        upper_border_lines.append(segment)
        lower_border_lines.append(segment)
        return Polygon(upper_border_lines, u_l, up, k=self.k), \
                Polygon(lower_border_lines, l_l, down, k = self.k)

    def score_split(self, segment, vertices):
        u_b_v, l_b_v = horizontal_split_vertices(vertices, segment)
        return abs(len(u_b_v) - len(l_b_v))

    def find_pretty_good_split(self, k = 3):
        """
        Checks all the pairs of vertices and finds the best pair of vertices with which to split
        this polygon with
        :return:
        """
        min_val = float("inf")
        max_segment = None

        vertices = self.get_border_vertices()
        pts = random.choices(vertices, k=k)
        for p1, p2 in itertools.combinations(pts):
            try:
                segment = Segment(to_line(p1, p2), p1[0], p2[0])
                tmp_val = self.score_split(segment, vertices)
                if min_val > tmp_val:
                    max_segment = segment
                    min_val = tmp_val
            except Exception:
                continue
        return max_segment


    def get_weight(self) -> float:
        return self.weight

    def to_complicated(self) -> bool:
        return False
        #return len(self.border_lines) > self.k


class PolyTree:

    def is_active(self, node):
        return node.get_weight() > self.min_weight

    def register_leaf(self, node):
        if self.is_active(node):
            self.active_set.add(node)

    def unregister(self, node):
        if node in self.active_set:
            self.active_set.remove(node)

    def get_heaviest(self):
        return max(self.get_leaves(), key=lambda x: x.pt_count())

    def __init__(self, weighted_lines,  points=list(), min_weight=-1, k = 100):
        li = -float("inf")
        ri = float("inf")
        self.root = Polygon(w_lines=[(Segment(l, li, ri), w) for l, w in weighted_lines], points=points, k = k)
        self.lines = set()
        self.min_weight = min_weight

        self.active_set = set()
        self.register_leaf(self.root)

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def get_leaves(self):
        stack = deque([self.root])

        all_traps = []
        while stack:
            curr_node = stack.pop()
            if curr_node.is_terminal():
                all_traps.append(curr_node)
            else:
                stack.append(curr_node.get_b())
                stack.append(curr_node.get_a())

        return all_traps

    def visualize_arrangement(self, ax, min_x, max_x, min_y, max_y):
        for poly in self.get_leaves():
            poly.visualize(ax, min_x, max_x, min_y, max_y)

    def visualize(self, file_name):
        def toOutputString(n):
            if isinstance(n, Segment_Node):
                name = "S"
            elif isinstance(n, Polygon):
                name = "T"
            return name

        def gen_label(n, l):
            if isinstance(n, Segment_Node):
                name = "Line(%f, %f, %f, %f, %d), %s" % (n.segment.a, n.segment.b, n.segment.xl,
                                                         n.segment.xr, n.segment.id, l)
            elif isinstance(n, Polygon):
                name = ("P(%d, %d")%(n.id, n.weight)
            return name
        stack = deque([(self.root, 0, "r")])
        graph = pydot.Dot(graph_type='graph', rankdir="LR")
        while stack:
            curr_node, nid, l = stack.pop()
            node = pydot.Node("%s_%d"% (toOutputString(curr_node), nid), label=gen_label(curr_node, l))
            graph.add_node(node)
            if not curr_node.is_terminal():
                if curr_node.get_r_or_a() is not None:
                    edge = pydot.Edge("%s_%d"% (toOutputString(curr_node), nid), "%s_%d" %
                                 (toOutputString(curr_node.get_r_or_a()), 2 * nid + 1))
                    graph.add_edge(edge)
                    stack.append((curr_node.get_r_or_a(), 2 * nid + 1, "r_a"))
                if curr_node.get_l_or_b() is not None:
                    edge = pydot.Edge("%s_%d" % (toOutputString(curr_node), nid), "%s_%d" %
                                  (toOutputString(curr_node.get_l_or_b()), 2 * nid + 2))

                    graph.add_edge(edge)
                    stack.append((curr_node.get_l_or_b(), 2 * nid + 2, "l_b"))
        graph.write_png(file_name + ".png")




    def insert_segment(self, parent_node, polygon, new_segment) -> Tuple[Polygon, Polygon]:
        upper, lower = polygon.horz_split(new_segment)
        new_s = Segment_Node(new_segment, upper, lower)
        if parent_node is None:
            self.root = new_s
        elif parent_node.get_b() == polygon:
            parent_node.set_b(new_s)
        else:
            parent_node.set_a(new_s)
        return upper, lower

    def insert_line(self, line, curr_node=None, merge=False):
        """
        Converts a line into many seperate segments.
        """
        for l in self.lines:
            if l.same_line(line):
                return
        else:
            self.lines.add(line)

        curr_node = self.root if curr_node is None else curr_node
        full_segment = Segment(line, -float("inf"), float("inf"))
        node_stack = deque([(curr_node, None, full_segment)])

        while node_stack:

            curr_node, parent_node, curr_segment = node_stack.pop()
            if curr_node.is_terminal():
                if self.is_active(curr_node):
                    upper, lower = self.insert_segment(parent_node, curr_node, curr_segment)
                    if upper.to_complicated():
                        node_stack.append((upper, curr_node, upper.find_pretty_good_split()))
                    if lower.to_complicated():
                        node_stack.append((lower, curr_node, lower.find_pretty_good_split()))

                    self.unregister(curr_node)
                    self.register_leaf(upper)
                    self.register_leaf(lower)
                else:
                    continue
            elif curr_node.crosses(curr_segment):
                upper_split, lower_split, upper_right = curr_segment.split(curr_node.segment)
                node_stack.append((curr_node.get_b(), curr_node, lower_split))
                node_stack.append((curr_node.get_a(), curr_node, upper_split))

            elif curr_node.segment_a(curr_segment):
                node_stack.append((curr_node.get_a(), curr_node, curr_segment))
            else:
                node_stack.append((curr_node.get_b(), curr_node, curr_segment))

    def zone(self, line, curr_node=None):
        """
        Returns all the cells this line crosses.
        """

        curr_node = self.root if curr_node is None else curr_node
        full_segment = Segment(line, -float("inf"), float("inf"))
        node_stack = deque([(curr_node, None, full_segment)])
        while node_stack:
            curr_node, parent_node, curr_segment = node_stack.pop()
            if curr_node.is_terminal():
                yield curr_node
            elif curr_node.segment.same_line(curr_segment):
                continue
            elif curr_node.crosses(curr_segment):
                upper_split, lower_split, upper_right = curr_segment.split(curr_node.segment)
                node_stack.append((curr_node.get_b(), curr_node, lower_split))
                node_stack.append((curr_node.get_a(), curr_node, upper_split))
            elif curr_node.segment_a(curr_segment):
                node_stack.append((curr_node.get_a(), curr_node, curr_segment))
            else:
                node_stack.append((curr_node.get_b(), curr_node, curr_segment))

def compute_cutting(test_set, weight_map, points, r):
    total_weight = sum(weight_map[l] for l in test_set)
    min_weight = total_weight / r
    # Cutting size
    #min_pt_count =  len(points) / (12 * r * r) - 1
    lines = weighted_shuffle(test_set, [weight_map[l] / float(total_weight) for l in test_set])
    #random.shuffle(test_set)
    tree = PolyTree([(l, weight_map[l]) for l in lines],
                       points,
                       min_weight=min_weight)
    for l in lines:
        if len(tree.active_set) == 0:
            return tree
        tree.insert_line(l)
    return tree





#visualize_cuttings3(1000, -4, 4, -4, 4, 2, 16)