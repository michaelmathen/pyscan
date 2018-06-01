
import pprint
import math

GLOBAL_TOL = 1e-08

def isclose(a, b):
    return abs(a - b) <= GLOBAL_TOL * max(abs(a), abs(b))


def approx_eq(a, b):
    """
    a < b
    :param a:
    :param b:
    :return:
    """
    return math.isclose(a, b)


def approx_above(a, b):
    """
    a < b
    :param a:
    :param b:
    :return:
    """
    if approx_eq(a, b):
        return False
    else:
        return a < b


def approx_below(a, b):
    """
    a > b
    :param a:
    :param b:
    :return:
    """
    if approx_eq(a, b):
        return False
    else:
        return a > b


def approx_eq_above(a, b):
    """
    a < b
    :param a:
    :param b:
    :return:
    """
    if approx_eq(a, b):
        return True
    else:
        return a <= b


def approx_eq_below(a, b):
    """
    a < b
    :param a:
    :param b:
    :return:
    """
    if approx_eq(a, b):
        return True
    else:
        return a >= b


def restrict(x, min_x, max_x):
    return min(max_x, max(x, min_x))

class Line:

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __getitem__(self, index):
        if index == 0:
            return self.a
        elif index == 1:
            return self.b
        else:
            raise IndexError("Lines only have 0 and 1 indices")


    def crossing_pt(self, line):
        return (self.x_intercept(line), self.y_intercept(line))

    def x_intercept(self, line):
        return (self.b - line.b) / (line.a - self.a)

    def y_intercept(self, line):
        return (line.b * self.a - self.b * line.a) / (self.a - line.a)

    def y_evaluate(self, y):
        return (y - self.b) / self.a

    def evaluate(self, x):
        if approx_eq(self.a, 0):
            """
            This is so if the line is used as the top of a trapezoid it will still be defined 
            even if inf * 0.
            """
            return self.b
        return x * self.a + self.b

    def eval_at_boundaries(self, line, xl, xr):
        yll = line.evaluate(xl)
        ylr = line.evaluate(xr)
        ysl = self.evaluate(xl)
        ysr = self.evaluate(xr)
        return ysl, ysr, yll, ylr

    def is_segment(self):
        return False

    def above_closed_interval(self, line, xl, xr):
        """
        Checks to see if the passed in line is below
        self over the interval. Can cross the line at the boundary.
        This will evaluate the crossing to be overlapping if they are very close.
        """
        if xl == -math.inf and xr == math.inf:
            return line.a == self.a and line.b <= self.b
        else:
            if xl == -math.inf:
                return approx_eq_above(line.evaluate(xr), self.evaluate(xr)) and approx_eq_below(line.a, self.a)
            elif xr == math.inf:
                return approx_eq_above(line.evaluate(xl), self.evaluate(xl)) and approx_eq_above(line.a, self.a)
            else:
                ysl, ysr, yll, ylr = self.eval_at_boundaries(line, xl, xr)
                return approx_eq_above(yll, ysl) and approx_eq_above(ylr, ysr)

    def below_closed_interval(self, line, xl, xr):
        return line.above_closed_interval(self, xl, xr)

    # def below_closed_interval(self, line, xl, xr):
    #     ysl, ysr, yll, ylr = self.eval_at_boundaries(line, xl, xr)
    #     if xl == -math.inf and xr == math.inf:
    #         if line.a == self.a:
    #             return approx_eq_below(line.b, self.b)
    #         else:
    #             return False
    #     elif xl == -math.inf:
    #         return approx_eq_below(ylr, ysr) and approx_eq_above(line.a, self.a)
    #     elif xr == math.inf:
    #         return approx_eq_below(yll, ysl) and approx_eq_below(line.a, self.a)
    #     else:
    #         return approx_eq_below(yll, ysl) and approx_eq_below(ylr, ysr)

    def above_interval(self, line, xl, xr):
        ysl, ysr, yll, ylr = self.eval_at_boundaries(line, xl, xr)
        if xl == -math.inf and xr == math.inf:
            if line.a == self.a:
                return approx_above(line.b, self.b)
            else:
                return False
        elif xl == -math.inf:
            return approx_above(ylr, ysr) and approx_eq_below(line.a, self.a)
        elif xr == math.inf:
            return approx_above(yll, ysl) and approx_eq_above(line.a, self.a)
        else:
            return approx_above(yll, ysl) and approx_above(ylr, ysr)

    def below_interval(self, line, xl, xr):
        return line.above_interval(self, xl, xr)

    def interval_crossed_by_closed(self, line, xl, xr):
        if approx_eq(line.a, self.a):
            return False
        x_val = self.x_intercept(line)
        return approx_eq_above(xl, x_val) and approx_eq_above(x_val, xr)

    def interval_crossed_by(self, line, xl, xr):
        if approx_eq(line.a, self.a):
            return False
        x_val = self.x_intercept(line)
        return approx_above(xl, x_val) and approx_above(x_val, xr)

    def same_line(self, other):
        return approx_eq(self.a, other.a) and approx_eq(self.b, other.b)

    def pt_eq_below(self, pt):
        y = self.evaluate(pt[0])
        return approx_eq_above(pt[1], y)

    def pt_eq_below_exact(self, pt):
        y = self.evaluate(pt[0])
        return pt[1] <= y

    def pt_eq_above(self, pt):
        y = self.evaluate(pt[0])
        return approx_eq_above(y, pt[1])

    def __eq__(self, other):
        if isinstance(other, Line):
            return other.a == self.a and other.b == self.b
        return False

    def to_dual(self):
        return (self.a, -self.b)

    def __hash__(self):
        return hash((self.a, self.b))

    def visualize(self, ax, min_x, max_x):
        ax.plot([min_x, max_x], [self.evaluate(min_x), self.evaluate(max_x)])


class Segment(Line):

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, line, x_1, x_2):
        super(Segment, self).__init__(line.a, line.b)
        self.xl = x_1
        self.xr = x_2

    def hsplit(self, x):
        """
        Splits the edge. Does not update the simplices though
        that the child edges point at. Need to update that later.
        """
        e1 = Segment(self, self.xl, x)
        e2 = Segment(self, x, self.xr)
        return e1, e2

    def is_segment(self):
        return True

    def split(self, line):
        """
        Splits the edge. Does not update the simplices though
        that the child edges point at. Need to update that later.
        """
        x_mid = self.x_intercept(line)
        if approx_eq(x_mid, self.xl) or approx_eq(x_mid, self.xr):
            if self.above_closed(line):
                return self, None, approx_eq(x_mid, self.xr)
            else:
                return None, self, approx_eq(x_mid, self.xr)

        e1 = Segment(self, self.xl, x_mid)
        e2 = Segment(self, x_mid, self.xr)
        if e1.above_closed(line):
            return e1, e2, True
        else:
            return e2, e1, False

    def simple_split(self, line):
        x_mid = self.x_intercept(line)
        if approx_eq(x_mid, self.xl) or approx_eq(x_mid, self.xr):
            if self.above_closed(line):
                return self, None
            else:
                return None, self
        e1 = Segment(self, self.xl, x_mid)
        e2 = Segment(self, x_mid, self.xr)
        if e1.above_closed(line):
            return e1, e2
        else:
            return e2, e1

    def below(self, line):
        return self.below_interval(line, self.xl, self.xr)

    def above(self, line):
        return self.above_interval(line, self.xl, self.xr)

    def above_closed(self, line):
        return self.above_closed_interval(line, self.xl, self.xr)

    def below_closed(self, line):
        return self.below_closed_interval(line, self.xl, self.xr)

    def crossed_by(self, line):
        return self.interval_crossed_by(line, self.xl, self.xr)

    def crossed_by_closed(self, line):
        return self.interval_crossed_by_closed(line, self.xl, self.xr)

    """
    Used to do segment vs. segment comparisons. Only tests on the overlapping portions 
    of the segment. Concludes false if there is no overlap.
    """
    def segment_overlap(self, seg):
        return self.xl < seg.xr and seg.xl < self.xr

    def overlap_region(self, seg):
        return (max(self.xl, seg.xl), min(self.xr, seg.xr))

    def below_segment(self, seg):
        xl, xr = self.overlap_region(seg)
        if xl >= xr:
            return False
        return self.below_interval(seg, xl, xr)

    def above_segment(self, line):
        xl, xr = self.overlap_region(line)
        if xl >= xr:
            return False
        return self.above_interval(line, xl, xr)

    def above_closed_segment(self, line):
        xl, xr = self.overlap_region(line)
        if xl >= xr:
            return False
        return self.above_closed_interval(line, xl, xr)

    def below_closed_segment(self, line):
        xl, xr = self.overlap_region(line)
        if xl >= xr:
            return False
        return self.below_closed_interval(line, xl, xr)

    def crossed_by_segment(self, line):
        xl, xr = self.overlap_region(line)
        if xl >= xr:
            return False
        return self.interval_crossed_by(line, xl, xr)

    def crossed_by_closed_segment(self, line):
        xl, xr = self.overlap_region(line)
        if approx_above(xr, xl):
            return False
        return self.interval_crossed_by_closed(line, xl, xr)

    @property
    def left_vertex(self):
        return self.xl, self.evaluate(self.xl)

    @property
    def right_vertex(self):
        return self.xr, self.evaluate(self.xr)

    def y_intersects(self, a, b):
        y1 = self.xl * a + b
        y2 = self.xr * a + b
        return y1, y2

    def contains_x_closed(self, x_val):
        return approx_eq_above(self.xl, x_val) and approx_eq_above(x_val, self.xr)

    def contains_x(self, x_val):
        return approx_above(self.xl, x_val) and approx_above(x_val, self.xr)

    @property
    def yl(self):
        return self.evaluate(self.xl)

    @property
    def yr(self):
        return self.evaluate(self.xr)

    def same_segment(self, other):
        if other is None:
            return False
        elif not isinstance(other, Edge):
            return False
        else:
            return approx_eq(other.a, self.a) and approx_eq(other.b, self.b) and approx_eq(other.xl, self.xl) and approx_eq(other.xr, self.xr)


    #This is useful for determining if this line segment is to the right side or left side of its intersection
    # note this doesn't make sense if xl and xr are -inf and inf
    def mostly_to_right(self, x):
        return (self.xl + self.xr) / 2.0 > x

    def mostly_to_left(self, x):
        return (self.xl + self.xr) / 2.0 < x

    def contained(self, pt):
        return self.pt_eq_above(pt) and self.pt_eq_below(pt)

    def to_left(self, x):
        return approx_above(self.xr, x)

    def to_right(self, x):
        return approx_above(x, self.xl)

    def to_left_closed(self, x):
        return approx_eq_above(self.xr, x)

    def to_right_closed(self, x):
        return approx_eq_above(x, self.xl)

    def pt_eq_below_seg(self, pt):
        if pt[0] > self.xr or pt[0] < self.xl:
            return False
        y = self.evaluate(pt[0])
        return approx_eq_above(pt[1], y)

    def pt_eq_above_seg(self, pt):
        if pt[0] > self.xr or pt[0] < self.xl:
            return False
        y = self.evaluate(pt[0])
        return approx_eq_above(y, pt[1])


def l_wedge(segment):
    up = None
    down = None
    x_mid = segment.a
    if math.isfinite(segment.xr):
        up = Segment(to_dual_line(segment.left_vertex), -math.inf, x_mid)
    if math.isfinite(segment.xl):
        down = Segment(to_dual_line(segment.right_vertex), -math.inf, x_mid)
    return Wedge(up, down)


def r_wedge(segment):
    up = None
    down = None
    x_mid = segment.a
    if math.isfinite(segment.xl):
        up = Segment(to_dual_line(segment.right_vertex), x_mid, math.inf)
    if math.isfinite(segment.xr):
        down = Segment(to_dual_line(segment.left_vertex), x_mid, math.inf)
    return Wedge(up, down)


class Wedge:

    def __init__(self, up_segment, down_segment):
        self.up_segment = up_segment
        self.down_segment = down_segment

    def contains(self, segment):
        if self.down_segment is not None and self.up_segment is not None:
            return (segment.above_closed(self.down_segment) and
                    segment.below_closed(self.up_segment))
        elif self.down_segment is not None:
            return segment.above_closed_segment(self.down_segment)
        elif self.up_segment is not None:
            return segment.below_closed_segment(self.up_segment)
        else:
            return True

    def crosses(self, segment):
        if self.down_segment is not None and self.up_segment is not None:
            return segment.crossed_by_closed_segment(self.up_segment) or \
                    segment.crossed_by_closed_segment(self.down_segment)
        elif self.down_segment is not None:
            return segment.crossed_by_closed_segment(self.down_segment)
        elif self.up_segment is not None:
            return segment.crossed_by_closed_segment(self.up_segment)
        else:
            return False

    def contains_pt(self, pt):
        if self.down_segment is not None and self.up_segment is not None:
            return self.up_segment.pt_eq_below_seg(pt) and self.down_segment.pt_eq_above_seg(pt)
        elif self.down_segment is not None:
            return self.down_segment.pt_eq_above_seg(pt)
        elif self.up_segment is not None:
            return self.up_segment.pt_eq_below_seg(pt)
        else:
            return True

    def above_closed(self, segment):
        """
        Must run contains first
        :param segment:
        :return:
        """
        if self.down_segment is not None:
            return self.down_segment.above_closed(segment)
        elif self.up_segment is not None:
            return self.up_segment.above_closed(segment)
        else:
            return False


    def split(self, segment):
        if self.up_segment.crossed_by_segment(segment):
            up_us, down_us = self.up_segment.simple_split(segment)
        else:
            up_us, down_us = (self.up_segment, self.up_segment)

        if self.down_segment.crossed_by_segment(segment):
            up_ds, down_ds = self.down_segment.simple_split(segment)
        else:
            up_ds, down_ds = (self.down_segment, self.up_segment)
        return Wedge(up_us, up_ds), Wedge(down_us, down_ds)


def to_line(p1, p2):
    return Line((p1[1] - p2[1]) / (p1[0] - p2[0]), p1[1] - (p1[1] - p2[1]) / (p1[0] - p2[0]) * p1[0])


def to_dual_line(pt):
    return Line(pt[0], -pt[1])


def to_dual_pt(line):
    return line.a, -line.b


def det2(a1, a2, b1, b2):
    return a1 * b2 - a2 * b1


def deduplicate_points(pts):
    if not pts:
        return []
    s_pts_x = sorted(pts, key=lambda x: x[0])

    curr_p = s_pts_x[0]
    groups = [[curr_p]]
    for p in s_pts_x[1:]:
        if approx_eq(curr_p[0], p[0]):
            groups[-1].append(p)
        else:
            curr_p = p
            groups.append([p])
    not_duplicate = []
    for g in groups:
        g.sort(key=lambda x: x[1])
        c_p = g[0]
        not_duplicate.append(c_p)
        for p in g[1:]:
            if not approx_eq(c_p[1], p[1]):
                not_duplicate.append(p)
                c_p = p

    return not_duplicate





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

def close_point(p1, p2):
    return approx_eq((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0.0)


def remove_close_points(edge_points):
    unique_pts = []
    for pt in edge_points:
        for pt_2 in unique_pts:
            if close_point(pt, pt_2):
               break
        else:
            unique_pts.append(pt)

    return unique_pts


def dis(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def triangle_area(p1, p2, p3):
    a = dis(p1, p2)
    b = dis(p2, p3)
    c = dis(p3, p1)
    s = 1 / 2.0 * (a + b + c)
    if s * (s - a) * (s - b) * (s - c) <= 0:
        return 0.0
    return math.sqrt(s * (s - a) * (s - b) * (s - c))


def simplex_center(boundary_points):
    x_avg = 0
    y_avg = 0
    for (x, y) in boundary_points:
        x_avg += x
        y_avg += y
    return x_avg / len(boundary_points), y_avg / len(boundary_points)


def weighted_avg(pts, ws):
    x_avg = 0
    y_avg = 0
    w_total = 0
    for (x, y), w in zip(pts, ws):
        x_avg += x * w
        y_avg += y * w
        w_total += w
    if w_total <= 0:
        return (0, 0)
    return (x_avg / w_total, y_avg / w_total)



def simplex_area_center(simplex, min_x, max_x, min_y, max_y):
    edge_points = []
    for e in simplex.edges:
        if e.xl > max_x or e.xr < min_x:
            continue
        x1 = restrict(e.xl, min_x, max_x)
        y1 = e.evaluate(x1)

        x2 = restrict(e.xr, min_x, max_x)
        y2 = e.evaluate(x2)

        if y1 > max_y < y2 or y1 < min_y > y2:
            continue

        y1 = restrict(y1, min_y, max_y)
        y2 = restrict(y2, min_y, max_y)
        edge_points.append((x1, y1))
        edge_points.append((x2, y2))
    edge_points = remove_close_points(edge_points)
    if len(edge_points) == 0:
        return (0, 0)
    elif len(edge_points) == 1:
        return edge_points[0]
    elif len(edge_points) == 2:
        return (edge_points[0][0] + edge_points[1][0]) / 2, (edge_points[0][1] + edge_points[1][1]) / 2
    triangle_point = edge_points[0]

    centers = []
    weights = []
    for e1, e2 in zip(edge_points[1:], edge_points[2:]):
        w = triangle_area(triangle_point, e1, e2)
        s_c = simplex_center([triangle_point, e1, e2])
        centers.append(s_c)
        weights.append(w)
    return weighted_avg(centers, weights)



