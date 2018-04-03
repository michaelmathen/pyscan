from Cuttings import approx_above, \
    approx_eq, Line
import pprint
import itertools
import random
from collections import deque
from SeidelTree import Segment, to_line
import math
import pydot
import numpy.random as npr

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


def horizontal_split_vertices(points, segment):
    up = []
    down = []
    for p in points:
        if segment.pt_eq_below(p):
            down.append(p)
        else:
            up.append(p)
    return up, down


def restrict(x, min_x, max_x):
    return min(max_x, max(x, min_x))


poly_id = 0

style = {"simplex_color": 'k',
         "simplex_alpha": .4,
         "simplex_line_thickness": 3,
         "edge_color": 'k',
         "line_color": 'k',
         "line_thickness": 1,
         "zone_line_color": "b",
         "zone_line_thickness": 4}


def visualize_edge(ax, e, min_x, max_x, min_y, max_y, c, linewidth):
    if e.xl < min_x > e.xr or e.xl > max_x < e.xr:
        return

    x1 = restrict(e.xl, min_x, max_x)
    y1 = e.evaluate(x1)

    x2 = restrict(e.xr, min_x, max_x)
    y2 = e.evaluate(x2)

    ax.plot([x1, x2], [y1, y2], c, linewidth=linewidth)


def split_lines(line, segments, ekey=lambda x: x, pkey=lambda x, y: x):
    u_b_l = []
    l_b_l = []
    for l_p in segments:
        l = ekey(l_p)
        if line.same_line(l):
            continue
        elif l.crossed_by(line):
            u_s, l_s = l.simple_split(line)
            u_b_l.append(pkey(u_s, l_p))
            l_b_l.append(pkey(l_s, l_p))
        elif l.above_closed(line):
            u_b_l.append(l_p)
        else:
            l_b_l.append(l_p)
    return u_b_l, l_b_l


class Polygon(Node):

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, border_lines=list(), w_lines=list(), points=list(), k=8):
        self.border_lines = border_lines
        self.w_lines = w_lines
        self.points = points
        self.weight = sum(w for _, w in w_lines)
        global poly_id
        self.id = poly_id
        poly_id += 1
        self.k = k

    def visualize(self, ax, min_x, max_x, min_y, max_y):
        if len(self.border_lines) <= 1:
            return
        pts = self.get_border_vertices()
        for e in self.border_lines:
            try:
                visualize_edge(ax, e, min_x, max_x, min_y, max_y, style["simplex_color"], style["simplex_line_thickness"])
            except ZeroDivisionError:
                continue
        local_pts = [p for p in pts if min_x <= p[0] <= max_x and min_y <= p[1] <= max_y]
        if local_pts:
            xs, ys = zip(*local_pts)
            ax.scatter(xs, ys)

        for l, _ in self.w_lines:
             visualize_edge(ax, l, min_x, max_x, min_y, max_y, style["line_color"], style["line_thickness"])

    def get_border_vertices(self):
        l_pts = [l.left_vertex for l in self.border_lines]
        r_pts = [l.right_vertex for l in self.border_lines]
        l_cycle = itertools.cycle(l_pts)
        next(l_cycle)

        border_pts = []
        for p1, p2 in zip(l_cycle, r_pts):
            if approx_eq(p1[0], p2[0]) and approx_eq(p1[1], p2[1]):
                border_pts.append(p1)
            else:
                border_pts.append(p1)
                border_pts.append(p2)

        return border_pts

    def get_vertices(self):
        return self.get_border_vertices()

    def get_points(self):
        return self.points

    def get_lines(self):
        return [l for l, _ in self.w_lines]

    def pt_count(self):
        return len(self.points)

    def is_terminal(self):
        return True

    def horz_split(self, segment):
        up, down = horizontal_split_vertices(self.points, segment)
        u_l, l_l = split_lines(segment, self.w_lines, ekey=lambda x: x[0], pkey=lambda x, y: (x, y[1]))
        u_b_l, l_b_l = split_lines(segment, self.border_lines)
        u_b_l.append(segment)
        l_b_l.append(segment)
        return Polygon(u_b_l, u_l, up, k=self.k), \
                Polygon(l_b_l, l_l, down, k=self.k)

    def score_split(self, segment, vertices):
        u_b_v, l_b_v = horizontal_split_vertices(vertices, segment)
        return abs(len(u_b_v) - len(l_b_v))

    def find_pretty_good_split_v(self):
        """
        Checks all the pairs of vertices and finds the best pair of vertices with which to split
        this polygon with
        :return:
        """
        min_val = float("inf")
        max_segment = None

        vertices = self.get_border_vertices()
        vertices = [p for p in vertices if not (math.isinf(p[0]) or math.isinf(p[1]))]
        p1 = random.choice(vertices)
        vertices.remove(p1)
        for p2 in vertices:
            try:
                segment = Segment(to_line(p1, p2), min(p1[0], p2[0]), max(p1[0], p2[0]))
                tmp_val = self.score_split(segment, vertices)
                if min_val > tmp_val:
                    max_segment = segment
                    min_val = tmp_val
            except ZeroDivisionError:
                pass
        return max_segment

    def find_pretty_good_split_l(self):
        vals = [w for _, w in self.w_lines]
        total_w = sum(vals)
        p = [w / total_w for w in vals]
        segments = npr.choice([l for l, _ in self.w_lines], p=p)
        return segments

    def get_weight(self) -> float:
        return self.weight

    def to_complicated(self) -> bool:
        #return False
        return len(self.border_lines) > self.k


class PolyTree:

    def is_active(self, node):
        return node.get_weight() > self.min_weight

    def get_heaviest(self):
        return max(self.get_leaves(), key=lambda x: x.pt_count())

    def __init__(self, weighted_lines,  points=list(), min_weight=-1, k = 8, seg_cutting=True):
        li = -float("inf")
        ri = float("inf")
        if seg_cutting:
            w_lines = []
            for l, w in weighted_lines:
                if l.is_segment():
                    w_lines.append((l, w))
                    #print("got here")
                else:
                    w_lines.append((Segment(l, li, ri), w))

            self.root = Polygon(w_lines=w_lines, points=points, k = k)
        else:

            self.root = Polygon(w_lines=[(Segment(l, li, ri), w) for l, w in weighted_lines], points=points, k = k)
        self.lines = set()
        self.min_weight = min_weight


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

    def insert_segment(self, parent_node, polygon, new_segment):
        upper, lower = polygon.horz_split(new_segment)
        new_s = Segment_Node(new_segment, upper, lower)
        if parent_node is None:
            self.root = new_s
        elif parent_node.get_b() == polygon:
            parent_node.set_b(new_s)
        else:
            parent_node.set_a(new_s)
        return upper, lower, new_s

    def cutting_greedy(self):
        """
        Converts a line into many seperate segments.
        """
        node_stack = deque([(self.root, None)])

        while node_stack:
            curr_node, parent_node = node_stack.pop()
            if self.is_active(curr_node) or curr_node.to_complicated():
                if curr_node.to_complicated():
                   segment = curr_node.find_pretty_good_split_v()
                else:
                    segment = curr_node.find_pretty_good_split_l()
                upper, lower, new_parent = self.insert_segment(parent_node, curr_node, segment)
                node_stack.append((upper, new_parent))
                node_stack.append((lower, new_parent))
            else:
                continue

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


def compute_cutting(test_set, weight_map, points, r, k=8):
    total_weight = sum(weight_map[l] for l in test_set)
    min_weight = total_weight / r
    # Cutting size
    tree = PolyTree([(l, weight_map[l]) for l in test_set],
                    points,
                    min_weight=min_weight, k=k)
    tree.cutting_greedy()
    return tree






# def random_box_line():
#     left_y = random.random()
#     right_y = random.random()
#     return Line(left_y - right_y, left_y)
#

# def test_set_points(pts, t):
#     test_set = []
#     rand_pts = random.sample(pts, int(math.sqrt(2 * t) + 1))
#     for p1, p2 in itertools.combinations(rand_pts, 2):
#         if len(test_set) >= t:
#             break
#         test_set.append(to_line(p1, p2))
#     return test_set
#
#
# def test_set_lines(pts, t):
#     test_set = []
#     rand_pts = random.sample(pts, 2 * t)
#
#     for p1, p2 in zip(rand_pts[:t], rand_pts[t:]):
#         test_set.append(to_line(p1, p2))
#     return test_set
# #
# #
# pts = [(random.random(), random.random()) for i in range(1000)]
# lines = test_set_lines(pts, 25)
# tree = compute_cutting_greedy(lines, {l: 1 for l in lines}, pts, 5, 5)
# matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
# f, ax = plt.subplots()
# tree.visualize_arrangement(ax, -1, 2, -1, 2)
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
# plt.show()



#visualize_cuttings3(1000, -4, 4, -4, 4, 2, 16)