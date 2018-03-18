from Cuttings import *
from collections import deque
import pydot
import numpy as np
import itertools
import time

id = 0
class Segment(Line):

    def __repr__(self):
        parameters = vars(self)
        non_recursive = {}
        for el in parameters:
            if isinstance(parameters[el], Edge) or isinstance(parameters[el], Simplex):
                non_recursive[el] = parameters[el].id
            else:
                non_recursive[el] = parameters[el]
        return type(self).__name__ + "(**" + pprint.pformat(non_recursive, indent=4, width=1) + ")"

    def __init__(self, line, x_1, x_2):
        super(Segment, self).__init__(line.a, line.b)
        self.xl = x_1
        self.xr = x_2
        global id
        self.id = id
        id += 1

    def hsplit(self, x):
        """
        Splits the edge. Does not update the simplices though
        that the child edges point at. Need to update that later.
        """
        e1 = Segment(self, self.xl, x)
        e2 = Segment(self, x, self.xr)
        return e1, e2

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



class Node:
    def is_terminal(self):
        return False

    def point_l_or_b(self, pt):
        return None

    def point_r_or_a(self, pt):
        return None

    def segment_l_or_b(self, segment):
        return None

    def segment_r_or_a(self, segment):
        return None

    def crosses(self, segment):
        return None

class Non_Terminal(Node):

    def __init__(self, left, right):
        self.right = right
        self.left = left

    def set_l_or_b(self, left):
        self.left = left

    def set_r_or_a(self, right):
        self.right = right

    def get_l_or_b(self):
        return self.left

    def get_r_or_a(self):
        return self.right

class Vertical_Node(Non_Terminal):
    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"
    def __init__(self, x, left=None, right=None):
        super(Vertical_Node, self).__init__(left, right)
        self.x = x



    def point_l_or_b(self, pt):
        return approx_eq_above(pt[0], self.x)

    def point_r_or_a(self, pt):
        return self.splitter.pt_is_right(self.x, pt[0])

    def hinted_r_or_a(self, seg, right):
        if right:
            if approx_eq(self.x, seg.right_vertex[0]):
                return False
            else:
                return approx_above(self.x, seg.right_vertex[0])
        else:
            if approx_eq(self.x, seg.left_vertex[0]):
                return True
            else:
                return approx_above(self.x, seg.left_vertex[0])

    def hinted_l_or_b(self, seg, right):
        return not self.hinted_r_or_a(seg, right)

    def segment_l_or_b(self, segment):
        """
        Determine if the segment is to the left of this vertical node.
        :param segment:
        :return:
        """
        return segment.to_left_closed(self.x)

    def segment_r_or_a(self, segment):
        return segment.to_right_closed(self.x)

    def crosses(self, seg):
        return approx_above(seg.xl, self.x) and \
            approx_above(self.x, seg.xr)


class Segment_Node(Non_Terminal):

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, segment, up=None, down=None):
        super(Segment_Node, self).__init__(down, up)
        self.segment = segment


    def point_l_or_b(self, pt):
        return self.segment.pt_eq_above(pt)

    def point_r_or_a(self, pt):
        return self.segment.pt_eq_below(pt)


    def hinted_r_or_a(self, seg, right):

        if right:
            if self.segment.contained(seg.right_vertex):
                return self.segment.below_closed_interval(seg, float("-inf"), seg.xr)
            else:
                return self.segment.pt_eq_above(seg.right_vertex)
        else:
            if self.segment.contained(seg.left_vertex):
                return self.segment.below_closed_interval(seg, seg.xl, float("inf"))
            else:
                return self.segment.pt_eq_above(seg.left_vertex)

    def hinted_l_or_b(self, seg, right):
        return not self.hinted_r_or_a(seg, right)

    def segment_l_or_b(self, seg):
        """
        Determine if the segment is below this node.
        :param segment:
        :return:
        """
        return seg.below_closed_interval(self.segment, seg.xl, seg.xr)

    def segment_r_or_a(self, seg):
        return seg.above_closed_interval(self.segment, seg.xl, seg.xr)

    def crosses(self, seg):
        """
        Has to intersect the interior of the segment.
        :param seg:
        :return:
        """
        x_v = seg.x_intercept(self.segment)
        return approx_above(self.segment.xl, x_v) and \
               approx_above(x_v, self.segment.xr) and \
               approx_above(seg.xl, x_v) and \
               approx_above(x_v, seg.xr)


DEBUG = True

trap_id = 0
class Trapezoid(Node):
    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, top_line=None, bottom_line=None, left_x=-float("inf"), right_x=float("inf"), lines=list(),
                 weights=list(), points=list()):
        self.top_line = Line(0, float("inf")) if top_line is None else top_line

        self.bottom_line = Line(0, float("-inf")) if bottom_line is None else bottom_line
        self.left_x = left_x
        self.right_x = right_x
        self.lines = []

        self.points = points
        # for p in points:
        #     #if self.inside_closed(p):
        #     self.points.append(p)

        self.weights = []
        for l, w in zip(lines, weights):
            if l.same_line(top_line) or l.same_line(bottom_line):
                continue
            if self.line_crosses(l):
                self.lines.append(l)
                self.weights.append(w)

        self.weight = sum(self.weights)
        global trap_id
        self.id = trap_id
        trap_id += 1

        global DEBUG
        if DEBUG:
            if approx_above(self.top_line.evaluate(self.left_x), self.bottom_line.evaluate(self.left_x)):
                print(self)
            elif approx_above(self.top_line.evaluate(self.right_x), self.bottom_line.evaluate(self.right_x)):
                print(self)

    def get_points(self):
        return self.points

    def ptCount(self):
        return len(self.points)

    def is_terminal(self):
        return True

    def visualize(self, ax, min_x, max_x, min_y, max_y, viz_lines=False, color='r'):

        if self.left_x > max_x or self.right_x < min_x:
            return

        if viz_lines:
            c = color
        else:
            c = next(ax._get_lines.prop_cycler)['color']
        print(self)
        print(self.top_line.evaluate(self.left_x), self.bottom_line.evaluate(self.left_x))
        print(self.top_line.evaluate(self.right_x), self.bottom_line.evaluate(self.right_x))

        if max_x >= self.left_x >= min_x:
            ax.vlines(self.left_x, min(max_y, self.top_line.evaluate(self.left_x)),
                        max(min_y, self.bottom_line.evaluate(self.left_x)), color="b")
        if max_x >= self.right_x >= min_x:
            ax.vlines(self.right_x, min(max_y, self.top_line.evaluate(self.right_x)),
                    max(min_y, self.bottom_line.evaluate(self.right_x)), color="b")
        mnx = max(self.left_x, min_x)
        mxx = min(self.right_x, max_x)
        ax.plot([mnx, mxx],
                [self.top_line.evaluate(mnx), self.top_line.evaluate(mxx)],
                linewidth=style["simplex_line_thickness"], color="b")
        ax.plot([mnx, mxx],
                [self.bottom_line.evaluate(mnx), self.bottom_line.evaluate(mxx)],
                linewidth=style["simplex_line_thickness"], color="b")
        if viz_lines:
            for l in self.lines:
                ax.plot([min_x, max_x], [l.evaluate(min_x), l.evaluate(max_x)], color=color)

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

        return Trapezoid(self.top_line, segment, self.left_x, self.right_x, self.lines, self.weights, up), \
                Trapezoid(segment, self.bottom_line, self.left_x , self.right_x, self.lines, self.weights, down)


    def vert_split(self, x):
        left = []
        right = []
        for p in self.points:
            if p[0] < x:
                left.append(p)
            else:
                right.append(p)
        return Trapezoid(self.top_line, self.bottom_line, self.left_x, x, self.lines, self.weights, left), \
               Trapezoid(self.top_line, self.bottom_line, x, self.right_x, self.lines, self.weights, right)


    def getWeight(self):
        return self.weight

    def inside_closed(self, pt):
        return approx_eq_above(self.left_x, pt[0]) and approx_eq_above(pt[0], self.right_x) and \
                self.top_line.pt_eq_below(pt) and self.bottom_line.pt_eq_above(pt)

    def crosses(self, segment):
        if self.inside_closed(segment.left_vertex) or self.inside_closed(segment.right_vertex):
            return True
        else:
            # check to see if it crosses the interior
            # TODO
            pass

    def line_crosses(self, line):
        # print(line.below_interval(self.top_line, self.left_x, self.right_x))
        # print(line.above_interval(self.bottom_line, self.left_x, self.right_x))
        # print(line.interval_crossed_by_closed(self.top_line, self.left_x, self.right_x))
        # print(line.interval_crossed_by_closed(self.bottom_line, self.left_x, self.right_x))

        if self.top_line.b == float("inf") and self.bottom_line.b == -float("inf"):
            return True
        elif self.top_line.b == float("inf"):

            return line.above_interval(self.bottom_line, self.left_x, self.right_x) or \
                   line.interval_crossed_by(self.bottom_line, self.left_x, self.right_x)
        elif self.bottom_line.b == float("-inf"):
            return line.below_interval(self.top_line, self.left_x, self.right_x) or \
                   line.interval_crossed_by(self.top_line, self.left_x, self.right_x)
        else:
            return (line.below_interval(self.top_line, self.left_x, self.right_x) and
                    line.above_interval(self.bottom_line, self.left_x, self.right_x)) or \
                    line.interval_crossed_by_closed(self.top_line, self.left_x, self.right_x) or \
                    line.interval_crossed_by_closed(self.bottom_line, self.left_x, self.right_x)




class Seidel_Tree:

    def is_active(self, node):
        return node.get_weight() > self.min_weight


    def add_act_heav(self, node):
        if self.is_active(node):
            self.active_set.add(node)
        self.heavy_set.add(node)

    def rem_act_heav(self, node):
        if node in self.active_set:
            self.active_set.remove(node)
        if node in self.heavy_set:
            self.heavy_set.remove(node)

    def get_heaviest(self):
        return max(self.heavy_set, key=lambda x: x.pt_count())

    def __init__(self, crossing_lines=list(), weights=None,  points=list(), min_weight=-1, min_p_count=-1):
        if weights is None:
            weights = [1 for i in range(len(crossing_lines))]
        self.root = Trapezoid(lines=crossing_lines, weights=weights, points=points)
        self.lines = set()
        self.min_weight = min_weight
        self.min_p_count = min_p_count

        self.heavy_set = set()
        self.active_set = set()
        self.add_act_heav(self.root)

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def visualize(self, file_name):
        def toOutputString(n):
            if isinstance(n, Segment_Node):
                name = "S"
            elif isinstance(n, Trapezoid):
                name = "T"
            elif isinstance(n, Vertical_Node):
                name = "V"
            return name

        def gen_label(n, l):
            if isinstance(n, Segment_Node):
                name = "Line(%f, %f, %f, %f, %d), %s" % (n.segment.a, n.segment.b, n.segment.xl,
                                                         n.segment.xr, n.segment.id, l)
            elif isinstance(n, Trapezoid):
                name = ("T(%d, %d, Line(%f, %f), \n " +
                       " Line(%f, %f)) \n" +
                        " lx=%f, rx=%f \n" +
                        "lty=%f, lby=%f, rty=%f rby=%f \n" +
                        "%s ")%(n.id, n.weight,
                                         n.top_line.a, n.top_line.b,
                                         n.bottom_line.a, n.bottom_line.b,
                                         n.left_x, n.right_x,
                                         n.top_line.evaluate(n.left_x), n.bottom_line.evaluate(n.left_x),
                                         n.top_line.evaluate(n.right_x),n.bottom_line.evaluate(n.right_x), l)
            elif isinstance(n, Vertical_Node):
                name = "V(%f), %s" % (n.x, l)
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

    def visualize_arrangement(self, ax, min_x, max_x, min_y, max_y):
        stack = deque([(self.root, [])])

        def max_val(segs, x_val):
            return max((s.evaluate(x_val) for s in segs))

        def min_val(segs, x_val):
            return min((s.evaluate(x_val) for s in segs))

        while stack:
            curr_node, segments = stack.pop()
            if not curr_node.is_terminal():
                if isinstance(curr_node, Segment_Node):
                    e = curr_node.segment
                    x1 = restrict(e.xl, min_x, max_x)
                    y1 = e.evaluate(x1)

                    x2 = restrict(e.xr, min_x, max_x)
                    y2 = e.evaluate(x2)
                    #m_y = min(min_y, y1, y2)
                    #mx_y = max(max_y, y1, y2)

                    ax.plot([x1, x2], [y1, y2], 'o')
                    ax.plot([x1, x2], [y1, y2])
                    segments += [e]
                elif isinstance(curr_node, Vertical_Node):
                    top_y = max_val(segments, curr_node.x)
                    bottom_y = min_val(segments, curr_node.x)
                    if min_x < curr_node.x < max_x:
                        ax.vlines(curr_node.x, max(bottom_y, min_y), min(top_y, max_y))
                if curr_node.get_r_or_a() is not None:
                    stack.append((curr_node.get_r_or_a(), segments))
                if curr_node.get_l_or_b() is not None:
                    stack.append((curr_node.get_l_or_b(), segments))

    def visualize_trapezoids(self, ax, min_x, max_x, min_y, max_y, color_list):
        """
        Picks a few trapezoids and colors them based on the color list. If 4 colors are passed in then
        four trapezoids will be colored. Uses the
        :param ax:
        :param min_x:
        :param max_x:
        :param min_y:
        :param max_y:
        :param color_list:
        :return:
        """
        stack = deque([self.root])

        all_traps = []
        while stack:
            curr_node = stack.pop()
            if curr_node.is_terminal():
                all_traps.append(curr_node)
            else:
                stack.append(curr_node.get_l_or_b())
                stack.append(curr_node.get_r_or_a())

        random.seed(1)
        random.shuffle(all_traps)

        # ix = 0
        # curr_node = all_traps[6]
        # print(len(all_traps))
        # curr_node.visualize(ax, min_x, max_x, min_y, max_y, viz_lines=False, color=color_list[0])

        for curr_node in all_traps[:len(color_list)]:
            curr_node.visualize(ax, min_x, max_x, min_y, max_y, viz_lines=True, color=color_list[ix])
        """
        for curr_node in all_traps[len(color_list):]:
            curr_node.visualize(ax, min_x, max_x, min_y, max_y)
        """

    def locate_pt(self, pt, curr_node=None):
        """
        Traces down the tree to find the location to insert the point.
        :param pt:
        :return:
        """
        curr_node = self.root if curr_node is None else curr_node
        parent_node = None
        while not curr_node.is_terminal():
            parent_node = curr_node
            if curr_node.point_l_or_b(pt):
                curr_node = curr_node.get_l_or_b()
            else:
                curr_node = curr_node.get_r_or_a()
        return curr_node, parent_node

    def locate_segment_end(self, line, right, curr_node=None):
        """
        Traces down the tree to find the location to insert the point. Uses the
        associated line as a hint to find the correct trapezoid.
        :param pt:
        :return:
        """
        curr_node = self.root if curr_node is None else curr_node
        parent_node = None
        while not curr_node.is_terminal():
            parent_node = curr_node
            if curr_node.hinted_r_or_a(line, right):
                curr_node = curr_node.get_r_or_a()
            else:
                curr_node = curr_node.get_l_or_b()
        return curr_node, parent_node

    def locate_segment(self, segment, curr_node=None):
        """
        generates tuples where the tuples consist of terminal nodes that we
        cross and their direct parents.
        Note this allows for segments that cross other segments.
        """

        curr_node = self.root if curr_node is None else curr_node
        node_stack = deque([(curr_node, None)])
        while node_stack:
            curr_node, parent_node = node_stack.pop()
            if curr_node.is_terminal():
                yield curr_node, parent_node
            elif curr_node.crosses(segment):
                node_stack.append((curr_node.get_r_or_a(), curr_node))
                node_stack.append((curr_node.get_l_or_b(), curr_node))
            elif curr_node.segment_r_or_a(segment):
                node_stack.append((curr_node.get_r_or_a(), curr_node))
            else:
                node_stack.append((curr_node.get_l_or_b(), curr_node))


    def insert_splitter(self, line, right, curr_node=None, parent_node=None):
        """
        :return:
        """

        curr_node = self.root if curr_node is None else curr_node

        if right:
            x_sp = line.xr
        else:
            x_sp = line.xl


        l_parent_node = None
        while not curr_node.is_terminal():
            l_parent_node = curr_node
            if isinstance(l_parent_node, Vertical_Node) and approx_eq(l_parent_node.x, x_sp):
                return
            elif curr_node.hinted_r_or_a(line, right):
                curr_node = curr_node.get_r_or_a()
            else:
                curr_node = curr_node.get_l_or_b()
        terminal_node = curr_node

        if not self.is_active(terminal_node):
            return

        if l_parent_node is not None:
            parent_node = l_parent_node

        if terminal_node.left_x > x_sp or terminal_node.right_x < x_sp:
            #self.visualize("Error_occured%f_%d"%(x_sp, terminal_node.id))
            #print("got here")
            return

        left_term, right_term = terminal_node.vert_split(x_sp)
        self.rem_act_heav(terminal_node)
        self.add_act_heav(left_term)
        self.add_act_heav(right_term)

        if parent_node is None:
            self.root = Vertical_Node(x_sp, left_term, right_term)
        elif parent_node.get_l_or_b() == terminal_node:
            parent_node.set_l_or_b(Vertical_Node(x_sp, left_term, right_term))
        else:
            parent_node.set_r_or_a(Vertical_Node(x_sp, left_term, right_term))

    def insert_line(self, line, curr_node=None):
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
                    upper, lower = curr_node.horz_split(curr_segment)
                    self.rem_act_heav(curr_node)
                    self.add_act_heav(upper)
                    self.add_act_heav(lower)

                    new_s = Segment_Node(curr_segment, upper, lower)
                    if parent_node is None:
                        self.root = new_s
                    elif parent_node.get_l_or_b() == curr_node:
                        parent_node.set_l_or_b(new_s)
                    else:
                        parent_node.set_r_or_a(new_s)
                else:
                    continue
            elif curr_node.crosses(curr_segment):
                if isinstance(curr_node, Segment_Node):
                    # upper right determines if we are inserting the right upper vertex or the left upper vertex
                    upper_split, lower_split, upper_right = curr_segment.split(curr_node.segment)
                    if upper_right:
                        self.insert_splitter(lower_split, False, curr_node=curr_node.get_l_or_b(),
                                             parent_node=curr_node)
                        self.insert_splitter(upper_split, True, curr_node=curr_node.get_r_or_a(),
                                             parent_node=curr_node)
                    else:
                        self.insert_splitter(lower_split, True, curr_node=curr_node.get_l_or_b(),
                                             parent_node=curr_node)
                        self.insert_splitter(upper_split, False, curr_node=curr_node.get_r_or_a(),
                                             parent_node=curr_node)

                    node_stack.append((curr_node.get_l_or_b(), curr_node, lower_split))
                    node_stack.append((curr_node.get_r_or_a(), curr_node, upper_split))
                elif isinstance(curr_node, Vertical_Node):
                    left_s, right_s = curr_segment.hsplit(curr_node.x)
                    node_stack.append((curr_node.get_l_or_b(), curr_node, left_s))
                    node_stack.append((curr_node.get_r_or_a(), curr_node, right_s))

            elif curr_node.segment_r_or_a(curr_segment):
                node_stack.append((curr_node.get_r_or_a(), curr_node, curr_segment))
            else:
                node_stack.append((curr_node.get_l_or_b(), curr_node, curr_segment))

    # def insert_segment(self, seg, curr_node=None, min_weight=-1):
    #     """
    #     The segment is assumed to not cross other line segments.
    #     :return:
    #     """
    #     curr_node = self.root if curr_node is None else curr_node
    #     for (n, p) in self.locate_segment(seg, curr_node=curr_node):
    #         if n.getWeight() <= min_weight:
    #             continue
    #         upper, lower = n.horz_split(seg)
    #         new_s = Segment_Node(seg, upper, lower)
    #         if p is None:
    #             self.root = new_s
    #         elif p.get_l_or_b() == n:
    #             p.set_l_or_b(new_s)
    #         else:
    #             p.set_r_or_a(new_s)
    #
    # def insert_full_segment(self, seg, curr_node=None, min_weight=-1):
    #     """
    #     The segment is assumed to not cross other line segments.
    #     :return:
    #     """
    #     curr_node = self.root if curr_node is None else curr_node
    #     if seg.xl != float("-inf"):
    #         self.insert_splitter(seg, False, curr_node, min_weight=min_weight)
    #     if seg.xr != float("inf"):
    #         self.insert_splitter(seg, True, curr_node, min_weight=min_weight)
    #     self.insert_segment(seg, curr_node, min_weight=min_weight)


def to_line(p1, p2):
    return Line((p1[1] - p2[1]) / (p1[0] - p2[0]), p1[1] - (p1[1] - p2[1]) / (p1[0] - p2[0]) * p1[0])


def compute_cutting(test_set, weight_map, points, r):
    total_weight = sum(weight_map[l] for l in test_set)
    min_weight = total_weight / r
    tree = Seidel_Tree(test_set, [weight_map[l] for l in test_set], points, min_weight=min_weight)
    for l in test_set:
        if len(tree.active_set) == 0:
            return tree
        tree.insert_line(l)
    return tree

def compute_slow_cutting(test_set, weight_map, points, r):
    total_weight = sum(weight_map[l] for l in test_set)
    min_weight = total_weight / r

    current_traps = [Trapezoid(lines=test_set, weights=[weight_map[l] for l in test_set], points=points)]
    for l in test_set:
        for t in current_traps:
            if t.line_crosses(l) and t.getWeight() > min_weight:
               pass


    #return tree

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


def partitions(pts, r, output_size, z=10):
    t = int(math.sqrt(r) + 1)
    #print(t)

    final_cells = []
    cell_queue = deque([pts])
    while cell_queue:

        curr_pts = cell_queue.pop()

        #rand_pts = random.sample(curr_pts, min(z, len(curr_pts)))
        test_set = []
        for i in range(z):
            p1, p2 = random.sample(curr_pts, 2)
            test_set.append(to_line(p1, p2))

        weight_map = {line: 1 for line in test_set}

        pts_not_in_cells = curr_pts[:]
        #print(len(pts_not_in_cells), len(curr_pts) / (2 * r))
        print("divide cell %d" % (len(curr_pts),))
        while len(pts_not_in_cells) > len(curr_pts) / (r**2) and len(pts_not_in_cells) > int(len(pts) / output_size):
            tree = compute_cutting(test_set, weight_map, pts_not_in_cells, r)
            print("Points left over %d" % (len(pts_not_in_cells),))
            #cells = {}
            if tree.heavy_set:
                cell = tree.get_heaviest()
                #print(cell)
                for l in cell.lines:
                    weight_map[l] *= 2
                lkup = set(cell.get_points())
                pts_not_in_cells = [p for p in pts_not_in_cells if p not in lkup]
                if len(cell.get_points()) < int(len(pts) / output_size) + 1:
                    final_cells.append(cell.get_points())
                else:
                    cell_queue.append(cell.get_points())
            else:
                random.shuffle(test_set)
                print("had to shuffle")
        if len(pts_not_in_cells) < int(len(pts) / output_size) + 1:
            final_cells.append(pts_not_in_cells)

    pts_Weights = [(random.sample(cell, 1)[0], len(cell) / float(len(pts))) for cell in final_cells]
    pts, weights = zip(*pts_Weights)
    return list(pts), list(weights)


def partition_tree(pts, r, output_size, z=10):
    t = int(math.sqrt(r) + 1)
    #print(t)

    tree_root = PartTree(Trapezoid())
    cell_queue = deque([(pts, tree_root)])

    while cell_queue:

        curr_pts, root = cell_queue.pop()

        rand_pts = random.sample(curr_pts, min(z, len(curr_pts)))
        test_set = []
        for i in range(z):
            p1, p2 = random.sample(curr_pts, 2)
            test_set.append(to_line(p1, p2))

        weight_map = {line: 1 for line in test_set}

        pts_not_in_cells = curr_pts[:]
        #print(len(pts_not_in_cells), len(curr_pts) / (2 * r))
        print("divide cell %d"%(len(curr_pts),))
        while len(pts_not_in_cells) > len(curr_pts) / (r**2) and len(pts_not_in_cells) > int(len(pts) / output_size) + 1:
            print("Points left over %d"%(len(pts_not_in_cells),))
            tree = compute_cutting(test_set, weight_map, pts_not_in_cells, r)

            #cells = {}
            if tree.heavy_set:
                cell = tree.get_heaviest()
                tree_leaf = PartTree(cell)
                root.insert_children(tree_leaf)
                for l in cell.lines:
                    weight_map[l] *= 2
                lkup = set(cell.get_points())
                pts_not_in_cells = [p for p in pts_not_in_cells if p not in lkup]
                cell_queue.append((cell.get_points(), tree_leaf))
            else:
                random.shuffle(test_set)
                print("had to shuffle")
        print("next level")
    return tree_root

def line_error_test(big_samp, small_samp, weights, pt_num):
    test_pool = big_samp + small_samp
    for i in range(pt_num):
        [p1, p2] = random.sample(test_pool, 2)
        line = to_line(p1, p2)

        r_1 = sum(line.pt_eq_below(p) for p in big_samp) / float(len(big_samp))

        r_2 = sum(weights[i] for i, p in enumerate(small_samp) if line.pt_eq_below(p)) / sum(weights[i] for i, p in enumerate(small_samp))

        yield abs(r_1 - r_2)


def error_plot(big_sample, l_s, h_s, k, r=4, z=100):

    errors = []
    sizes = []
    for i in np.linspace(l_s, h_s, k):
        s_size = int(i + .5)
        sampled_pts, sampled_weights = partitions(big_sample, r, s_size, z=z)
        #print(sampled_weights)
        max_error = 0
        for e in line_error_test(big_sample, sampled_pts, sampled_weights, 100):
            max_error = max(max_error, e)
        errors.append(max_error)
        sizes.append(s_size)

    f, a = plt.subplots()
    a.scatter(sizes, errors)
    plt.show()


def plot_cutting_size(pts, l_s, h_s, k, r=4):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = []
        rand_pts = random.sample(pts, int(math.sqrt(i)))

        for p1, p2 in itertools.combinations(rand_pts, 2):
            test_set.append(to_line(p1, p2))
        weight_map = {t: 1 for t in test_set}
        tree = compute_cutting(test_set, weight_map, pts, r)
        trap_count.append(len(tree.heavy_set) / float(r * r))
        line_count.append(len(test_set))

    f, a = plt.subplots()
    a.scatter(line_count, trap_count)
    plt.show()

def plot_cutting_size2(pts, l_s, h_s, k, r=4):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = [random_box_line() for j in range(int(i))]
        weight_map = {t: 1 for t in test_set}
        tree = compute_cutting(test_set, weight_map, pts, r)
        trap_count.append(len(tree.heavy_set) / float(r * r))
        line_count.append(len(test_set))

    f, a = plt.subplots()
    a.scatter(line_count, trap_count)
    plt.show()



# matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
pts = [(random.random(), random.random()) for i in range(10000)]
#plot_cutting_size(pts, 16, 1000, 80)
#plot_cutting_size2(pts, 16, 1000, 80)
# #compute_cutting()
error_plot(pts, 10, 400, 40, r=2, z=50)


# tree = partition_tree(pts, 2, 20, z=1000)
# f, a = plt.subplots()
# # #
# tree.visualize(a, -5, 5, -5, 5, viz_lines=False, color = 'r')
# xs, ys = zip(*pts)
# a.scatter(xs, ys)
# plt.show()

# Spatial data download (Get the lidar data to use)


# big_sample = [(random.random(), random.random()) for i in range(1000)]
# times = []
# r_vals = []
# for i in range(1, 20):
#     s = time.time()
#     sampled_pts, sampled_weights = partitions(big_sample, 4, 200, z=(i * 2))
#     times.append(time.time() - s)
#     r_vals.append(2 * i)
#     print(2 * i)
#
# f, a = plt.subplots()
# a.scatter(r_vals, times)
# plt.show()

"""
big_sample = [(random.random(), random.random()) for i in range(1000)]
times = []
r_vals = []
for i in range(1, 16):
    s = time.time()
    sampled_pts, sampled_weights = partitions(big_sample, 2 * int(i), 200)
    times.append(time.time() - s)
    r_vals.append(2 * int(i))
    print(2 * i)

f, a = plt.subplots()
a.scatter(r_vals, times)
plt.show()
"""
# if __name__ == "__main__":
#     #arrange = Arrangement([random_box_line() for i in range(10)])
#     #assume the domain is a 1 by 1 square.
#     #construct the lines so they extend from the left side to the right side.
#     max_x = 2
#     min_x = -2
#     # tree = Seidel_Tree([random_box_line() for i in range(100)], [1 for i in range(100)], min_weight=30)
#     matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     pts = [(random.random(), random.random()) for i in range(1000)]
#     sampled_pts, sampled_weights = partitions(pts, 16, 400)
#
#     f, a = plt.subplots()
#     x_vals, y_vals = zip(*sampled_pts)
#     a.scatter(x_vals, y_vals)
#     plt.show()

# if __name__ == "__main__":
#     #arrange = Arrangement([random_box_line() for i in range(10)])
#     #assume the domain is a 1 by 1 square.
#     #construct the lines so they extend from the left side to the right side.
#     max_x = 2
#     min_x = -2
#     # tree = Seidel_Tree([random_box_line() for i in range(100)], [1 for i in range(100)], min_weight=30)
#     matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     # for i in range(50):
#     #     #f, (a1, a2) = plt.subplots(2, 1)
#     #     new_line = random_box_line()
#     #     print(new_line)
#     #     #arrange.insert(new_line)
#     #     tree.insert_line(new_line)
#     #     #tree.visualize("testing_%d"%(i,))
#     #     #print(tree)
#     lines = [random_box_line() for i in range(1000)]
#     weights = {l: 1 for l in lines}
#
#     tree = compute_cutting(lines, weights, 4)
#
#     f, a = plt.subplots()
#     #tree.visualize_trapezoids(a, -10, 10, -10, 10, color_list=['r'])
#     tree.visualize_arrangement(a, -5, 5, -5, 5)
#     plt.show()

# if __name__ == "__main__":
#     l = Line(**{   'a': -0.4976264579860905, 'b': 0.3526024541774099})
#     t = Trapezoid(top_line=Line(**{   'a': 0, 'b': float("inf")}),
#               bottom_line=Line(**{'a': -0.36399907743655024, 'b': 0.2901120445807751}),
#               left_x=-float("inf"), right_x=-0.8119035028412716, lines=[l], weights=[1])
#     f, a = plt.subplots()
#     t.visualize(a, -5, 5, -5, 5, viz_lines=True, color = 'r')
#     plt.show()
#     print(t.line_crosses(l))


# if __name__ == "__main__":
#     # arrange = Arrangement([random_box_line() for i in range(10)])
#     # assume the domain is a 1 by 1 square.
#     # construct the lines so they extend from the left side to the right side.
#     max_x = 2
#     min_x = -2
#     tree = Seidel_Tree([random_box_line() for i in range(10)], [1 for i in range(10)])
#     matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     lines = [Line(**{ 'a': 0.5360426951034875,
#                         'b': 0.8326691657946446}),
#                         Line(**{   'a': -0.0480690624799861,
#                             'b': 0.49243669185492167}),
#                         Line(**{   'a': 0.4680545566595835,
#                             'b': 0.8112910112858986}),
#                         Line(**{   'a': 0.3898850932855328,
#                             'b': 0.8150769095233331})]
#
#     i = 0
#     for new_line in lines[:-1]:
#         # f, (a1, a2) = plt.subplots(2, 1)
#         #new_line = random_box_line()
#         #print(new_line)
#         # arrange.insert(new_line)
#         tree.insert_line(new_line, min_weight=-1)
#         tree.visualize("testing_%d"%(i,))
#         i+=1
#         # print(tree)
#
#     f, a = plt.subplots()
#     tree.insert_line(lines[-1], min_weight=-1)
#
#     #tree.visualize_trapezoids(a, -5, 5, -5, 5, color_list=['r'])
#     tree.visualize_arrangement(a, -5, 5, -5, 5)
#     plt.show()

# if __name__ == "__main__":
#     #arrange = Arrangement([random_box_line() for i in range(10)])
#     #assume the domain is a 1 by 1 square.
#     #construct the lines so they extend from the left side to the right side.
#     max_x = 2
#     min_x = -2
#     tree = Seidel_Tree([random_box_line() for i in range(100)], [1 for i in range(100)])
#     matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     lines = [Line(**{'a': -0.09945156480785433,
#             'b': 0.5699475893302439}),
#         Line(**{'a': 0.20333947989923118,
#             'b': 0.41909911052666216}),
#         Line(**{'a': -0.37812378855716655,
#             'b': 0.2988252079696353})]
#     i = 0
#
#     for new_line in lines[:-1]:
#         #f, (a1, a2) = plt.subplots(2, 1)
#         #arrange.insert(new_line)
#         tree.insert_line(new_line, min_weight=10)
#         tree.visualize("testing_%d"%(i,))
#         i+=1
#         f, a = plt.subplots()
#         tree.visualize_arrangement(a, -5, 5, -5, 5)
#         plt.show()
#         #print(tree)
#
#     new_line = lines[-1]
#     tree.insert_line(new_line, min_weight=10)
#     tree.visualize("testing_%d"%(i,))
#     i+=1
#     f, a = plt.subplots()
#     tree.visualize_arrangement(a, -5, 5, -5, 5)
#     plt.show()
#
#     plt.show()