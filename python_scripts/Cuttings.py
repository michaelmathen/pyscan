import random
import matplotlib.pyplot as plt
import matplotlib
import pprint
import pickle
import math

ix = 0
GLOBAL_TOL = 1e-08

def weighted_shuffle(items, weights):
    rweights = [-random.random() ** (1.0 / w) for w in weights]
    order = sorted(range(len(items)), key=lambda i: rweights[i])
    return [items[i] for i in order]

def approx_eq(a, b):
    """
    a < b
    :param a:
    :param b:
    :return:
    """
    return math.isclose(a, b, rel_tol=GLOBAL_TOL)

def approx_above(a, b):
    """
    a < b
    :param a:
    :param b:
    :return:
    """
    if math.isclose(a, b, rel_tol=GLOBAL_TOL):
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
    if math.isclose(a, b, rel_tol=GLOBAL_TOL):
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
    if math.isclose(a, b, rel_tol=GLOBAL_TOL):
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
    if math.isclose(a, b, rel_tol=GLOBAL_TOL):
        return True
    else:
        return a >= b

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
        if xl == -float("inf") and xr == float("inf"):
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
        ysl, ysr, yll, ylr = self.eval_at_boundaries(line, xl, xr)
        if xl == -math.inf and xr == math.inf:
            if line.a == self.a:
                return approx_eq_below(line.b, self.b)
            else:
                return False
        elif xl == -math.inf:
            return approx_eq_below(ylr, ysr) and approx_eq_above(line.a, self.a)
        elif xr == math.inf:
            return approx_eq_below(yll, ysl) and approx_eq_below(line.a, self.a)
        else:
            return approx_eq_below(yll, ysl) and approx_eq_below(ylr, ysr)

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
        ysl, ysr, yll, ylr = self.eval_at_boundaries(line, xl, xr)
        if xl == -math.inf and xr == math.inf:
            if line.a == self.a:
                return approx_below(line.b, self.b)
            else:
                return False
        elif xl == -math.inf:
            return approx_below(ylr, ysr) and approx_eq_above(line.a, self.a)
        elif xr == math.inf:
            return approx_below(yll, ysl) and approx_eq_below(line.a, self.a)
        else:
            return approx_below(yll, ysl) and approx_below(ylr, ysr)

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

    def pt_eq_above(self, pt):
        y = self.evaluate(pt[0])
        return approx_eq_above(y, pt[1])

    def __eq__(self, other):
        if isinstance(other, Line):
            return other.a == self.a and other.b == self.b
        return False

    def __hash__(self):
        return hash((self.a, self.b))



ix_edge = 0


class Edge(Line):

    def __repr__(self):
        parameters = vars(self)
        non_recursive = {}
        for el in parameters:
            if isinstance(parameters[el], Edge) or isinstance(parameters[el], Simplex):
                non_recursive[el] = parameters[el].id
            else:
                non_recursive[el] = parameters[el]
        return type(self).__name__ + "(**" + pprint.pformat(non_recursive, indent=4, width=1) + ")"

    def __init__(self, line, x_1, x_2, upper_s=None, lower_s=None):
        super(Edge, self).__init__(line.a, line.b)
        self.xl = x_1
        self.xr = x_2
        self.right_edge = None
        self.left_edge = None
        self.up_simplex = upper_s
        self.bottom_simplex = lower_s
        global ix_edge
        self.id = ix_edge
        ix_edge += 1

    def split(self, line):
        """
        Splits the edge. Does not update the simplices though
        that the child edges point at. Need to update that later.
        """
        x_mid = self.x_intercept(line)
        if approx_eq(x_mid, self.xl) or approx_eq(x_mid, self.xr):
            if self.above_closed(line):
                return self, None
            else:
                return None, self

        e1 = Edge(self, self.xl, x_mid, None, None)
        e2 = Edge(self, x_mid, self.xr, None, None)
        if self.left_edge is not None:
            self.left_edge.right_edge = e1
        e1.left_edge = self.left_edge
        e1.right_edge = e2
        e2.left_edge = e1
        e2.right_edge = self.right_edge
        if self.right_edge is not None:
            self.right_edge.left_edge = e2

        e1.up_simplex = self.up_simplex
        e1.bottom_simplex = self.bottom_simplex
        e2.up_simplex = self.up_simplex
        e2.bottom_simplex = self.bottom_simplex

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

    def update_simplex(self, old_s, new_s):
        if self.up_simplex == old_s:
            self.up_simplex = new_s
        elif self.bottom_simplex == old_s:
            self.bottom_simplex = new_s

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

    def same_edge(self, other):
        if other is None:
            return False
        elif not isinstance(other, Edge):
            return False
        else:
            return other.a == self.a and other.b == self.b and other.xl == self.xl and other.xr == self.xr


class Splitter:
    def __init__(self, term_l, old_line, new_line):
        # vertex
        self.term_l = term_l
        self.old_line = old_line
        self.new_line = new_line

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def x_value(self):
        return self.old_line.x_intercept(self.new_line)

    def y_cross_value(self):
        """
        The y value corresponding to the crossing point of old_line and new_line
        """
        return self.old_line.y_intercept(self.new_line)

    def base_above(self, line):
        return approx_above(line.evaluate(self.x_value()), self.y_cross_value())

    def base_above_closed(self, line):
        return approx_eq_above(line.evaluate(self.x_value()), self.y_cross_value())

    def intersects(self, line):
        by = self.bottom_y()
        ty = self.top_y()
        my = line.evaluate(self.x_value())
        return approx_above(by, my) and approx_above(my, ty)

    def intersects_closed(self, line):
        c_y = line.evaluate(self.x_value())
        return approx_eq_below(self.top_y(), c_y) and approx_eq_below(c_y, self.bottom_y())

    def term_y(self):
        return self.term_l.evaluate(self.x_value())

    def intersect_y(self):
        return self.old_line.evaluate(self.x_value())

    def top_y(self):
        return max(self.term_y(), self.intersect_y())

    def bottom_y(self):
        return min(self.term_y(), self.intersect_y())

    def below_closed(self, line):
        """
        Is the splitter below the line
        """
        return approx_eq_above(self.top_y(), line.evaluate(self.x_value()))

    def above_closed(self, line):
        """
        Is the splitter above the line
        """
        return approx_eq_below(self.bottom_y(), line.evaluate(self.x_value()))

    def below(self, line):
        """
        Is the splitter below the line
        """
        return approx_above(self.top_y(), line.evaluate(self.x_value()))

    def above(self, line):
        """
        Is the splitter above the line
        """
        return approx_below(self.bottom_y(), line.evaluate(self.x_value()))

    def get_old_line(self):
        return self.old_line

    def get_new_line(self):
        return self.new_line

    def get_term_line(self):
        return self.term_l

    def cut_left(self):
        """
        Does the new line cut the old line to the left of the vertical line
        :return:
        """
        if self.term_y() < self.intersect_y():
            return self.old_line.a < self.new_line.a
        else:
            return self.new_line.a <= self.old_line.a

    def cut_right(self):
        return not self.cut_left()

    def pt_is_left(self, pt):
        return approx_eq_above(pt[0], self.x_value())

    def pt_is_right(self, pt):
        return approx_eq_above(self.x_value(), pt[0])



class Simplex:
    def __init__(self, crossing_lines=[]):
        self.edges = []
        self.splitters = []
        self.crossing_lines = [crossing_lines]
        global ix
        self.id = ix
        ix += 1



    def inside_point(self):
        if len(edges) > 1:
            distinct_edges = []
            for e in self.edges:
                if len(distinct_edges) >= 3:
                    return simplex_center(distinct_edges)
                distinct_edges.append(e.right_vertex)
                distinct_edges.append(e.left_vertex)
                remove_close_points(distinct_edges)
        else:
            return None

    # def inside(self, x, y):
    #     """
    #     simple slow inside implementation.
    #     :param x:
    #     :param y:
    #     :return:
    #     """
    #     crossing_pts = []
    #     l_test = Line(0, y)
    #     for e in self.edges:
    #         if e.crossed_by_closed(l_test):
    #             crossing_pts.append(e.crossed_pt(l_test))
    #     remove_close_points(crossing_pts)
    #
    #     if len(crossing_pts)

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"


    def print(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def get_new_edge(self, line, entrance_edge):
        """
        Makes new edge for the exit values. Note this does not set the upper and lower
        simplex.
        """
        exit_edge = None
        for edge in self.edges:
            if edge.crossed_by_closed(line) and not edge.same_edge(entrance_edge):
                exit_edge = edge

        if exit_edge is None:
            right_x = float("inf")
        else:
            right_x = line.x_intercept(exit_edge)

        if entrance_edge is None:
            left_x = -float("inf")
        else:
            left_x = line.x_intercept(entrance_edge)

        return Edge(line, left_x, right_x)

    def find_term_line(self, line, edge):
        v_x = line.x_intercept(edge)
        for e in self.edges:
            if e.contains_x(v_x) and not edge.same_line(e):
                return e
        return None

    def exit_edge(self, line, entrance_edge):
        for e in self.edges:
            if e.crossed_by_closed(line) and not e.same_edge(entrance_edge):
                return e
        return None

    def next_simplex(self, line, entrance_edge):
        exit_edge = self.exit_edge(line, entrance_edge)
        if exit_edge is None:
            return None
        else:
            if exit_edge.crossed_by_closed(line):
                if exit_edge.up_simplex == self:
                    return exit_edge, exit_edge.bottom_simplex
                else:
                    return exit_edge, exit_edge.up_simplex
            else:  # crosses the vertex
                if exit_edge.right_edge is None:
                    return None
                if exit_edge.up_simplex == self:
                    return exit_edge, exit_edge.right_edge.bottom_simplex
                else:
                    return exit_edge, exit_edge.right_edge.up_simplex

    def get_left(self):
        if not self.edges:
            return -float("inf")
        return min(e.xl for e in self.edges)

    def get_right(self):
        if not self.edges:
            return float("inf")
        return max(e.xr for e in self.edges)

    def insert_splitter(self, line, entrance_edge, exit_edge):

        # First insert the possibly two new splitters.
        # We do this by finding the possibly two new vertices created.
        # Create a splitter for each of these.
        # These will be the endpoints of the new edge.

        new_splitters = self.splitters[:]
        new_crossing_lines = self.crossing_lines[:]
        for edge in [entrance_edge, exit_edge]:
            if edge is not None:
                # could be the entrance or exit edge.
                term = self.find_term_line(line, edge)
                if term is None:  # the line extends up or down forever.
                    if edge == exit_edge:
                        if line.a > edge.a:
                            term = Line(0, -float('inf'))
                        elif line.a < edge.a:
                            term = Line(0, float('inf'))
                        else:
                            continue
                    else:
                        if line.a < edge.a:
                            term = Line(0, -float('inf'))
                        elif line.a > edge.a:
                            term = Line(0, float('inf'))
                        else:
                            continue
                new_splitter = Splitter(term, edge, line)

                for i, el in enumerate(new_splitters):
                    if new_splitter.x_value() < el.x_value():
                        new_splitters.insert(i, new_splitter)
                        break
                else:
                    i = len(new_crossing_lines) - 1
                    new_splitters.append(new_splitter)
                left = []
                right = []
                for l in new_crossing_lines[i]:
                    if new_splitter.intersects(l):
                        left.append(l)
                        right.append(l)
                    elif abs(new_splitter.term_y()) == float("inf"):
                        old = new_splitter.get_old_line()
                        i_x = old.x_intercept(l)
                        splitter_x = new_splitter.x_value()
                        if i_x < splitter_x:
                            left.append(l)
                        else:
                            right.append(l)
                    else:
                        old = new_splitter.get_old_line()
                        term = new_splitter.get_term_line()

                        wedge_point = old.x_intercept(term)
                        splitter_x = new_splitter.x_value()
                        intercepts = [l.x_intercept(term), l.x_intercept(old)]

                        if wedge_point < splitter_x:
                            intercepts = [x for x in intercepts if x > wedge_point]
                        else:
                            intercepts = [x for x in intercepts if x < wedge_point]
                        if len(intercepts) == 0:
                            print(wedge_point, splitter_x, [l.x_intercept(term), l.x_intercept(old)])
                        elif intercepts[0] < splitter_x:
                            left.append(l)
                        else:
                            right.append(l)
                new_crossing_lines[i] = right
                new_crossing_lines.insert(i, left)

        return new_splitters, new_crossing_lines

    def cut_splitters(self, line, new_splitters):

        #determine if the splitter is above or below the new line.
        m_upper_splitters = []
        m_lower_splitters = []
        # | /-\
        for splitter in new_splitters:
            if splitter.base_above(line) or splitter.above_closed(line):
                m_upper_splitters.append(splitter)
            else:
                m_lower_splitters.append(splitter)

        for splitter in new_splitters:
            if splitter.intersects(line):
                splitter.term_l = line

        return m_upper_splitters, m_lower_splitters

    def split_crossing_lines(self, line, sp, cross_lines):
        # Partiton the crossing lines by checking if they cross the new splitter
        # if they don't then they must reside only on 1 side or the other.
        # Check to see if they intersect the upper or lower lines in the right or left splitter.
        partition_points = [s.x_value() for s in sp] + [self.get_right()]
        upper_lines = []
        lower_lines = []
        if len(sp) == 0:
            return cross_lines, cross_lines

        cl = cross_lines[0]; lx = self.get_left(); rx = partition_points[0]
        if (sp[0].get_new_line().same_line(line) and sp[0].cut_left()) or sp[0].intersects(line):
                upper_lines.append([l for l in cl if l.above_interval(line, lx, rx) or
                                    line.interval_crossed_by(l, lx, rx)])
                lower_lines.append([l for l in cl if l.below_interval(line, lx, rx) or
                                    line.interval_crossed_by(l, lx, rx)])
        elif sp[0].above_closed(line):
            upper_lines.append(cl)
            lower_lines.append([])
        else:
            upper_lines.append([])
            lower_lines.append(cl)


        for lx, cl, s, rx in zip(partition_points[:-1], cross_lines[1:], sp, partition_points[1:]):
            if (s.get_new_line().same_line(line) and s.cut_right()) or s.intersects(line):
                upper_lines.append([l for l in cl if l.above_interval(line, lx, rx) or
                                    line.interval_crossed_by(l, lx, rx)])
                lower_lines.append([l for l in cl if l.below_interval(line, lx, rx) or
                                    line.interval_crossed_by(l, lx, rx)])
            elif s.above_closed(line):
                upper_lines.append(cl)
                lower_lines.append([])
            else:
                upper_lines.append([])
                lower_lines.append(cl)


        return upper_lines, lower_lines

    def merge_cells(self, line, sp, upper_lines, lower_lines):

        prev_upper_list = upper_lines[0] if upper_lines else []
        prev_lower_list = lower_lines[0] if lower_lines else []
        m_upper_lines = []
        m_lower_lines = []

        for splitter, u_line_cell, l_line_cell in zip(sp, upper_lines[1:], lower_lines[1:]):
            if splitter.base_above(line) or splitter.above_closed(line):
                prev_lower_list.extend(l_line_cell)
                m_upper_lines.append(list(set(prev_upper_list)))
                prev_upper_list = u_line_cell
            else:
                prev_upper_list.extend(u_line_cell)
                m_lower_lines.append(list(set(prev_lower_list)))
                prev_lower_list = l_line_cell
        m_upper_lines.append(list(set(prev_upper_list)))
        m_lower_lines.append(list(set(prev_lower_list)))
        return m_upper_lines, m_lower_lines

    def split_simplex(self, line, entrance_edge, u_entrance = None, l_entrance = None):
        upper_edges = []
        lower_edges = []

        exit_edge = self.exit_edge(line, entrance_edge)
        new_edge = self.get_new_edge(line, entrance_edge)
        if exit_edge is not None:
            u_exit, l_exit = exit_edge.split(line)
        else:
            u_exit = None; l_exit = None
        new_edge_flag = True
        if not self.edges:
            upper_edges.append(new_edge)
            lower_edges.append(new_edge)
        for edge in self.edges:
            if edge == entrance_edge:
                if l_entrance is not None:
                    lower_edges.append(l_entrance)
                if new_edge_flag:
                    upper_edges.append(new_edge)
                    lower_edges.append(new_edge)
                    new_edge_flag = False
                if u_entrance is not None:
                    upper_edges.append(u_entrance)
            elif edge == exit_edge:
                if u_exit is not None:
                    upper_edges.append(u_exit)
                if new_edge_flag:
                    upper_edges.append(new_edge)
                    lower_edges.append(new_edge)
                    new_edge_flag = False
                if l_exit is not None:
                    lower_edges.append(l_exit)
            elif edge.above_closed(line):
                upper_edges.append(edge)
            elif edge.below_closed(line):
                lower_edges.append(edge)
        upper_simplex = Simplex()
        lower_simplex = Simplex()
        new_edge.up_simplex = upper_simplex
        new_edge.bottom_simplex = lower_simplex
        for e in upper_edges:
            e.update_simplex(self, upper_simplex)
        for e in lower_edges:
            e.update_simplex(self, lower_simplex)
        upper_simplex.edges = upper_edges
        lower_simplex.edges = lower_edges

        """
        Makes a pass through the splitters. First inserts the possibly two
        new splitters. Secondly prunes lines that only fall above
        or below a division between two splitters. Finally
        merges splitters that are terminated by lines outside
        of the current simplex.
        """
        new_splitters, new_crossing_lines = self.insert_splitter(line, entrance_edge, exit_edge)
        upper_lines, lower_lines = self.split_crossing_lines(line, new_splitters, new_crossing_lines)
        print(len(upper_lines))
        upper_lines, lower_lines = self.merge_cells(line, new_splitters, upper_lines, lower_lines)
        m_upper_splitters, m_lower_splitters = self.cut_splitters(line, new_splitters)
        print(len(upper_lines), len(m_upper_splitters))
        print()
        upper_simplex.splitters = m_upper_splitters
        lower_simplex.splitters = m_lower_splitters
        upper_simplex.crossing_lines = upper_lines
        lower_simplex.crossing_lines = lower_lines

        return new_edge, u_exit, l_exit, upper_simplex, lower_simplex


class Trapezoid:

    def __init__(self):
        """
        A trapezoid consists of upper cells and lower cells.
        """
        self.left_splitter = None
        self.top_edge = None
        self.right_splitter = None
        self.bottom_edge = None
        self.left_cells = []
        self.right_cells = []
        self.top_cells = []
        self.bottom_cells = []







class Arrangement:
    def __init__(self, initial_lines=[]):
        self.left_edges = set()
        self.simplices = {Simplex(initial_lines)}

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def check_arrangement(self):
        all_edges = set()
        for s in self.simplices:
            for e in s.edges:
                if e.up_simplex is not None and e.up_simplex not in self.simplices:
                    print("simplex=%d, edge=%d, up_bad=%d"%(s.id, e.id, e.up_simplex.id))
                if e.bottom_simplex is not None and e.bottom_simplex not in self.simplices:
                    print("simplex=%d, edge=%d, up_bad=%d" % (s.id, e.id, e.bottom_simplex.id))
                all_edges.add(e)
        for s in self.simplices:
            for e in s.edges:
                if e.right_edge is not None and e.right_edge not in all_edges:
                    print("edge=%d, pointing to=%d"%(e.id, e.right_edge.id))


    def zone(self, line):
        if not self.left_edges:
            el = self.simplices.pop()
            self.simplices.add(el)
            yield None, el
        else:
            exit_edge = min(self.left_edges, key=lambda x: x.x_intercept(line))
            while not exit_edge.crossed_by_closed(line):
                exit_edge = exit_edge.right_edge
                if exit_edge is None:
                    print("Shouldn't be here")
                    return
            if line.a < exit_edge.a:
                curr_simplex = exit_edge.up_simplex
            else:
                curr_simplex = exit_edge.bottom_simplex
            yield None, curr_simplex
            entrance_edge = None
            while True:
                possibly_next = curr_simplex.next_simplex(line, entrance_edge)
                # print(possibly_next)
                if possibly_next is None:
                    return
                entrance_edge, curr_simplex = possibly_next
                yield possibly_next

    def insert(self, line):
        first_edge = None
        prev_edge = None
        uex = None
        lex = None
        zone = [el for el in self.zone(line)]
        for entrance_edge, simplex in zone:
            # print("Entrance " + str(entrance_edge))
            # print("Exit " + str(simplex.print()))
            # print()
            # print()
            u_e = uex
            l_e = lex
            new_edge, uex, lex, up_s, l_s = simplex.split_simplex(line, entrance_edge, u_e, l_e)
            self.simplices.remove(simplex)
            self.simplices.add(up_s)
            self.simplices.add(l_s)
            # plot_edge(new_edge, -1, 2, ax)
            # plot_edge(u_e, -1, 2, ax)
            # plot_edge(l_e, -1, 2, ax)
            # plot_edge(uex, -1, 2, ax)
            # plot_edge(lex, -1, 2, ax)
            # plt.show()
            # make sure that the left edges are kept current.
            if entrance_edge in self.left_edges:
                self.left_edges.remove(entrance_edge)
                if u_e is not None and l_e is not None:
                    if u_e.right_edge == l_e:
                        self.left_edges.add(u_e)
                    else:
                        self.left_edges.add(l_e)
                elif u_e is not None:
                    self.left_edges.add(u_e)
                elif l_e is not None:
                    self.left_edges.add(l_e)

            if first_edge is None:
                first_edge = new_edge
                prev_edge = first_edge
            else:
                prev_edge.right_edge = new_edge
                if new_edge is not None:
                    new_edge.left_edge = prev_edge
                prev_edge = new_edge
        self.left_edges.add(first_edge)


style = {"simplex_color": 'k',
         "simplex_alpha": .4,
         "simplex_line_thickness": 4,
         "edge_color": 'k',
         "line_color": 'k',
         "line_thickness": 1,
         "zone_line_color": "b",
         "zone_line_thickness": 4}


def plot_edge(edge, min_x, max_x, ax):
    """
    Plots the edge object
    :param edge:
    :param min_x:
    :param max_x:
    :param ax:
    :return:
    """
    if edge is None:
        return

    lx = max(edge.xl, min_x)
    rx = min(edge.xr, max_x)
    if lx > max_x or rx < min_x:
        return
    ax.plot([lx, rx], [edge.evaluate(lx), edge.evaluate(rx)], style["edge_color"])


def plot_line(line, min_x, max_x, ax):
    """
    Plots the line object
    :param line:
    :param min_x:
    :param max_x:
    :param ax:
    :return:
    """
    ax.plot([min_x, max_x], [line.evaluate(min_x), line.evaluate(max_x)], style["line_color"])


def plot_zone_line(line, min_x, max_x, ax):
    """
    Plots the line object
    :param line:
    :param min_x:
    :param max_x:
    :param ax:
    :return:
    """
    ax.plot([min_x, max_x], [line.evaluate(min_x), line.evaluate(max_x)], style["zone_line_color"],
            linewidth=style["zone_line_thickness"])


def restrict(x, min_x, max_x):
    return min(max_x, max(x, min_x))

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

def visualize_simplex(simplex, min_x, max_x, ax, ymin=-5, ymax=5, plot_lines = True):
    """
    Takes a simplex and walks the edge list. Draws each edge to a matplotlib drawing.
    :param simplex:
    :return:
    """
    m_y = float("inf")
    mx_y = float("-inf")
    any_in_range = False
    for e in simplex.edges:
        if e.xl < min_x > e.xr or e.xr > max_x < e.xl:
            continue
        any_in_range = True
        x1 = restrict(e.xl, min_x, max_x)
        y1 = e.evaluate(x1)

        x2 = restrict(e.xr, min_x, max_x)
        y2 = e.evaluate(x2)
        m_y = min(m_y, y1, y2)
        mx_y = max(mx_y, y1, y2)

        ax.plot([x1, x2], [y1, y2], style["simplex_color"], linewidth=style["simplex_line_thickness"])
    if not any_in_range:
        return
    x_c, y_c = simplex_area_center(simplex, min_x, max_x, m_y, mx_y)
    if x_c != min_x and x_c != max_x:
        ax.plot(x_c, y_c, 'o' )
        ax.text(x_c, y_c, str(simplex.id))

    for split in simplex.splitters:
        if max_x < split.x_value() or split.x_value() < min_x:
            continue
        if split.intersect_y() < m_y > split.term_y() or split.intersect_y() > mx_y < split.term_y():
            continue
        ty = split.term_y()
        iy = split.intersect_y()

        #ax.plot([split.x_value(), split.x_value()], [ty, iy])
        ax.vlines(split.x_value(), max(split.bottom_y(), ymin), min(split.top_y(), ymax))
    if plot_lines:
        for cell in simplex.crossing_lines:
            c = (random.random(), random.random(), random.random())
            for l in cell:
                ax.plot([min_x, max_x], [l.evaluate(min_x), l.evaluate(max_x)], c=c)
    #ax.plot((max(x_avg) + min(x_avg)) / 2, (max(y_avg) - min(y_avg)) / 2, 'x')#, str(simplex.id))
    # ax.fill(x_vals, y_vals, style["simplex_color"], alpha=style["simplex_alpha"])


def visualize_zone(arr, line, min_x, max_x, ax):
    for e, simplex in arr.zone(line):
        visualize_simplex(simplex, min_x, max_x, ax)
    plot_zone_line(line, min_x, max_x, ax)


def visualize_arrangement(arr, min_x, max_x, ax):
    visualize_lines(arr.left_edges, min_x, max_x, ax)

def visualize_lines(lines, min_x, max_x, ax):
    for l in lines:
        ax.plot([min_x, max_x], [l.evaluate(min_x), l.evaluate(max_x)], style["line_color"])


def random_box_line():
    left_y = random.random()
    right_y = random.random()
    return Line(left_y - right_y, left_y)


def visualize_arrangement_s(arr, min_x, max_x, ax, plot_lines=False):
    for s in arr.simplices:
        if s is not None:
            visualize_simplex(s, min_x, max_x, ax, plot_lines=plot_lines)
    label_arrangement(arr, min_x, max_x, ax)

def label_edge(edge, min_x, max_x, ax):
    x1 = restrict(edge.xl, min_x, max_x)
    y1 = edge.evaluate(x1)

    x2 = restrict(edge.xr, min_x, max_x)
    y2 = edge.evaluate(x2)

    ax.text((x1 + x2) / 2.0 , (y1 + y2) / 2.0, str(edge.id))

def label_simplices(s, min_x, max_x, ax):
    for e in s.edges:
        label_edge(e, min_x, max_x, ax)

def label_arrangement(arr, min_x, max_x, ax):
    edges = set()
    for s in arr.simplices:
        for e in s.edges:
            edges.add(e)

    for e in edges:
        label_edge(e, min_x, max_x, ax)




# if __name__ == "__main__":
#     max_x = 2
#     min_x = -1
#     f, (a1, a2) = plt.subplots(2, 1)
#     arrange = Arrangement()
#     l1 = Line(-0.6657358683439032, 0.029299023845995475)
#     l2 = Line(-0.3608683090015692, 0.22122812561026362)
#     arrange.insert(l1)
#     arrange.insert(l2)
#     visualize_arrangement(arrange, min_x, max_x, a1)
#     visualize_arrangement_s(arrange, min_x, max_x, a2)
#     plt.show()

# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     f, (a1, a2) = plt.subplots(2, 1)
#     arrange = Arrangement()
#     l1 = Line(0.002105977438810358, 0.5635184908196953)
#     l2 = Line(-0.36517199495130637, 0.1248506586629391)
#     l3 = Line(-0.8372357009616058, 0.07308605057553241)
#     arrange.insert(l1)
#     arrange.insert(l2)
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     #visualize_zone(arrange, l3, min_x, max_x, a2)
#     arrange.insert(l3)
#     visualize_arrangement(arrange, min_x, max_x, a1)
#     visualize_arrangement_s(arrange, min_x, max_x, a2)
#     a2.set_ylim([-3, 3])
#     plt.show()


# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     f, (a1, a2) = plt.subplots(2, 1)
#     arrange = Arrangement()
#     l1 = Line(0.8146451615701115, 0.97995435114778)
#     l2 = Line(-0.19332038537275786,0.07952872043352255)
#     l3 = Line(0.24059657208737162, 0.8530601610061491)
#     arrange.insert(l1)
#     arrange.insert(l2)
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     # visualize_zone(arrange, l3, min_x, max_x, a2)
#     arrange.insert(l3)
#     visualize_arrangement(arrange, min_x, max_x, a1)
#     visualize_arrangement_s(arrange, min_x, max_x, a2)
#     plt.show()



# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     f, (a1, a2) = plt.subplots(2, 1)
#     arrange = Arrangement()
#     lines = [Line(0.14884468689120656,0.6284582894908913),
#             Line(-0.7584834213850338, 0.2089412541839617),
#             Line(0.2943056873920782, 0.775972784417263)]
#             #Line(-0.338931373907027, 0.1842907936695246)]
#     for line in lines:
#         arrange.insert(line)
#     l_end = Line(-0.39914928816595274, 0.48188861882939105)
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     # visualize_zone(arrange, l3, min_x, max_x, a2)
#     #arrange.insert(l_end)
#     visualize_arrangement(arrange, min_x, max_x, a1)
#     visualize_arrangement_s(arrange, min_x, max_x, a2)
#     plt.show()

#def check_simplex_structure(arr):


# if __name__ == "__main__":
#     arrange = Arrangement([random_box_line() for i in range(20)])
#
#
#     #assume the domain is a 1 by 1 square.
#     #construct the lines so they extend from the left side to the right side.
#     max_x = 20
#     min_x = -20
#
#     matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     for i in range(10):
#         f, (a1, a2) = plt.subplots(2, 1)
#         new_line = random_box_line()
#         print(new_line)
#
#         #print(arrange)
#
#         #visualize_zone(arrange, new_line, min_x, max_x, a2)
#         arrange.insert(new_line)
#         visualize_arrangement_s(arrange, min_x, max_x, a1)
#
#
#         #visualize_arrangement_s(arrange, min_x, max_x, a2)
#         s = next(iter(arrange.simplices))
#         visualize_simplex(s, min_x, max_x, a2, plot_lines=True)
#         #print(s)
#         print()
#         print()
#         # for s in arrange.simplices:
#         #     print(s)
#         #     print()
#         # print()
#         # print()
#         # print()
#         for a_i in [a1, a2]:
#             a_i.set_xlim([min_x, max_x])
#             a_i.set_ylim([-3, 3])
#         plt.show()




"""

Line(**{   'a': -0.4893957882265083,
    'b': 0.20105278522340475})
Line(**{   'a': 0.11553119529761358,
    'b': 0.7844323075156047})
Line(**{   'a': 0.037134640852142464,
    'b': 0.6306839774421847})
Line(**{   'a': 0.31057577989759744,
    'b': 0.7433464169861587})

"""
# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     f, a1 = plt.subplots()
#     arrange = Arrangement()
#     lines = [Line(**{'a': -0.28002249578781524,
#                      'b': 0.04670219515585938})
#              , Line(**{'a': 0.06789914790966867,
#                      'b': 0.6770218081727364})
#              , Line(**{'a': -0.08642455382757375,
#                      'b': 0.29438253103340406})
#              , Line(**{'a': -0.003420918812250373,
#                      'b': 0.6343825278628497})]
#     for line in lines[:-1]:
#         arrange.insert(line)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#     label_arrangement(arrange, min_x, max_x, a1)
#     #l_end = Line(**{'a': -0.5824457065566974, 'b': 0.11903916420577654})
#     #l_end = Line(**{'a': -0.512094221280181, 'b': 0.0986132496986516})
#     #visualize_arrangement_s(arrange, min_x, max_x, a1)
#     l_end = lines[-1]
#     #visualize_zone(arrange, l_end, min_x, max_x, a1)
#     arrange.insert(l_end)
#     visualize_arrangement_s(arrange, min_x, max_x, a1)
#     plt.show()
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     # visualize_zone(arrange, l3, min_x, max_x, a2)
#     #arrange.insert(l_end)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#
#     #visualize_arrangement_s(arrange, min_x, max_x, a2)
#     #plt.show()


# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     f, a1 = plt.subplots()
#     arrange = Arrangement()
#     lines = [Line(0.6903139680905406, 0.7766000614458256),
#              Line(0.2471127832365867, 0.8804570217517127)]
#     for line in lines:
#         arrange.insert(line)
#     #visualize_arrangement_s(arrange, min_x, max_x, a1)
#     #label_arrangement(arrange, min_x, max_x, a1)
#
#     #l_end = Line(0.04264071968744698, 0.07750688115657534)
#     #visualize_arrangement_s(arrange, min_x, max_x, a1)
#     l_end = Line(-0.6626494625168489, 0.10553545514561435)
#     visualize_zone(arrange, l_end, min_x, max_x, a1)
#     plt.show()
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     # visualize_zone(arrange, l3, min_x, max_x, a2)
#     #arrange.insert(l_end)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#
#     #visualize_arrangement_s(arrange, min_x, max_x, a2)
#     #plt.show()

# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     matplotlib.rcParams['figure.figsize'] = [25.0, 25.0]
#
#     f, a1 = plt.subplots()
#     arrange = Arrangement()
#     lines = [Line(**{   'a': -0.06204898621984323,
#                         'b': 0.6085487637547395}),
#                         Line(**{   'a': -0.26184363232298846,
#                         'b': 0.4364408226629509}),
#                         Line(**{   'a': -0.8243131772394836,
#                         'b': 0.14888100302027707}),
#                         Line(**{   'a': 0.0648113588732635,
#                         'b': 0.8331150225923524}),
#                         Line(**{   'a': 0.4554075891321202,
#                         'b': 0.9592317992316349})]
#     for line in lines[:-1]:
#         arrange.insert(line)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#     label_arrangement(arrange, min_x, max_x, a1)
#     #l_end = Line(**{'a': -0.5824457065566974, 'b': 0.11903916420577654})
#     #l_end = Line(**{'a': -0.512094221280181, 'b': 0.0986132496986516})
#     #visualize_arrangement_s(arrange, min_x, max_x, a1)
#     l_end = lines[-1]
#     #visualize_zone(arrange, l_end, min_x, max_x, a1)
#     arrange.insert(l_end)
#     arrange.check_arrangement()
#     visualize_arrangement_s(arrange, min_x, max_x, a1)
#
#     plt.show()
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     # visualize_zone(arrange, l3, min_x, max_x, a2)
#     #arrange.insert(l_end)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#
#     #visualize_arrangement_s(arrange, min_x, max_x, a2)
#     #plt.show()



# if __name__ == "__main__":
#     max_x = 3
#     min_x = -3
#     matplotlib.rcParams['figure.figsize'] = [25.0, 25.0]
#
#     f, a1 = plt.subplots()
#     arrange = Arrangement()
#     lines = [Line(**{   'a': 0.5652948481586917,
#                 'b': 0.931025117704599}),
#             Line(**{   'a': -0.516417481230533,
#                 'b': 0.007590827422918167}),
#             Line(**{   'a': 0.2576243808977082,
#                 'b': 0.7766541977878596}),
#             Line(**{   'a': 0.4520475715133212,
#                 'b': 0.9229445088479115}),
#             Line(**{   'a': 0.4705553803152347,
#                 'b': 0.536457755244084})]
#     for line in lines[:-1]:
#         arrange.insert(line)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#
#     #l_end = Line(**{'a': -0.5824457065566974, 'b': 0.11903916420577654})
#     #l_end = Line(**{'a': -0.512094221280181, 'b': 0.0986132496986516})
#     #visualize_arrangement_s(arrange, min_x, max_x, a1)
#     l_end = lines[-1]
#     #visualize_zone(arrange, l_end, min_x, max_x, a1)
#     arrange.insert(l_end)
#     arrange.check_arrangement()
#     visualize_arrangement_s(arrange, min_x, max_x, a1)
#     label_arrangement(arrange, min_x, max_x, a1)
#     for s in arrange.simplices:
#         if s.id == 29:
#             print(s)
#     plt.show()
#     # f = open("test_file.pickle", 'w')
#     # pickle.dump(arrange, f)
#     # visualize_zone(arrange, l3, min_x, max_x, a2)
#     #arrange.insert(l_end)
#     #visualize_arrangement(arrange, min_x, max_x, a1)
#
#     #visualize_arrangement_s(arrange, min_x, max_x, a2)
#     #plt.show()

# if __name__ == "__main__":
#     max_x = 5
#     min_x = -5
#     matplotlib.rcParams['figure.figsize'] = [25.0, 25.0]
#
#     arrange = Arrangement([random_box_line() for i in range(8)])
#     lines = [Line(**{   'a': 0.06716988590351047,
#                         'b': 0.17757762985107906}),
#              Line(**{'a': -0.45931686931229443,
#                      'b': 0.37462296327674716}),
#              Line(**{'a': 0.39398499723011104,
#                      'b': 0.9811855638001208})]
#              # Line(**{'a': -0.4917764577689462,
#              #         'b': 0.38281185958504793})]
#
#     bef = 2
#     for line in lines[:bef]:
#         arrange.insert(line)
#
#     for line in lines[bef:]:
#         f, ((a1, a2), (a3, a4)) = plt.subplots(2, 2)
#         visualize_arrangement_s(arrange, min_x, max_x, a3)
#         plot_zone_line(line, min_x, max_x, a3)
#         for s in arrange.simplices:
#             if s.id == 3:
#                 visualize_simplex(s, min_x, max_x, a4, plot_lines=True)
#         arrange.insert(line)
#         for s in arrange.simplices:
#             if s.id == 11:
#
#                 visualize_simplex(s, min_x, max_x, a2, plot_lines=True)
#             if s.id == 12:
#                 visualize_simplex(s, min_x, max_x, a1, plot_lines=True)
#         plt.show()
#
#         print()


if __name__ == "__main__":
    arrange = Arrangement([random_box_line() for i in range(10)])


    #assume the domain is a 1 by 1 square.
    #construct the lines so they extend from the left side to the right side.
    max_x = 2
    min_x = -2

    matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
    for i in range(5):
        f, (a1, a2) = plt.subplots(2, 1)
        new_line = random_box_line()
        print(new_line)
        arrange.insert(new_line)

    for s in arrange.simplices:
        f, (a1, a2) = plt.subplots(2, 1)

        for a_i in [a1, a2]:
            a_i.set_xlim([min_x, max_x])
            a_i.set_ylim([-2, 2])
        visualize_arrangement_s(arrange, min_x, max_x, a1, plot_lines=False)

        visualize_simplex(s, min_x, max_x, a2, plot_lines=True)
        print(s)
        plt.savefig("fig_" + str(s.id) + ".png")