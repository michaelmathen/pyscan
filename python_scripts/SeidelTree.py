from Cuttings import *
from collections import deque
import pydot
import numpy as np
import itertools

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
        if approx_eq(seg.a, self.segment.a):
            return False
        x_v = seg.x_intercept(self.segment)
        return approx_above(self.segment.xl, x_v) and \
               approx_above(x_v, self.segment.xr) and \
               approx_above(seg.xl, x_v) and \
               approx_above(x_v, seg.xr)




DEBUG = True
trigered = False

trap_id = 0
class Trapezoid(Node):
    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, top_line=None, bottom_line=None, left_x=-float("inf"), right_x=float("inf"), w_lines=[], points=list()):
        self.top_line = Line(0, float("inf")) if top_line is None else top_line

        self.bottom_line = Line(0, float("-inf")) if bottom_line is None else bottom_line
        self.left_x = left_x
        self.right_x = right_x
        self.w_lines = []

        self.points = points
        # for p in points:
        #     #if self.inside_closed(p):
        #     self.points.append(p)
        self.weight = 0

        for l, w in w_lines:
            if l.same_line(self.top_line) or l.same_line(self.bottom_line):
                continue
            if self.line_crosses(l):
                self.w_lines.append((l, w))
                self.weight += w

        global trap_id
        self.id = trap_id
        trap_id += 1

        global DEBUG
        global trigered
        if DEBUG:
            if not approx_eq_above(self.bottom_line.evaluate(self.left_x),
                                                self.top_line.evaluate(self.left_x)):
                print(self)
                print("PROBLEM")
                trigered = True
            if not approx_eq_above(self.bottom_line.evaluate(self.right_x),
                                    self.top_line.evaluate(self.right_x)):
                print(self)
                print("PROBLEM")
                trigered = True

    def get_points(self):
        return self.points

    def get_lines(self):
        return [l for l, _ in self.w_lines]


    def pt_count(self):
        return len(self.points)

    def is_terminal(self):
        return True

    def visualize(self, ax, min_x, max_x, min_y, max_y, viz_lines=True, viz_points=False):

        if self.left_x > max_x or self.right_x < min_x:
            return

        if viz_lines:
            for l, _ in self.w_lines:
                ax.plot([min_x, max_x], [l.evaluate(min_x), l.evaluate(max_x)], color="k",
                        linewidth=style["line_thickness"])

        c = next(ax._get_lines.prop_cycler)['color']

        if max_x >= self.left_x >= min_x:
            ax.vlines(self.left_x, min(max_y, self.top_line.evaluate(self.left_x)),
                        max(min_y, self.bottom_line.evaluate(self.left_x)), color=c, linewidth=style["simplex_line_thickness"])
        if max_x >= self.right_x >= min_x:
            ax.vlines(self.right_x, min(max_y, self.top_line.evaluate(self.right_x)),
                    max(min_y, self.bottom_line.evaluate(self.right_x)), color=c, linewidth=style["simplex_line_thickness"])
        mnx = max(self.left_x, min_x)
        mxx = min(self.right_x, max_x)
        ax.plot([mnx, mxx],
                [self.top_line.evaluate(mnx), self.top_line.evaluate(mxx)],
                linewidth=style["simplex_line_thickness"], color=c)
        ax.plot([mnx, mxx],
                [self.bottom_line.evaluate(mnx), self.bottom_line.evaluate(mxx)],
                linewidth=style["simplex_line_thickness"], color=c)
        if viz_points:
            x_vals = []
            y_vals = []
            for x, y in self.points:
                if min_x <= x <= max_x and min_y <= y <= max_y:
                    x_vals.append(x)
                    y_vals.append(y)

            ax.scatter(x_vals, y_vals)


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

        # up_w_lines = []
        # down_w_lines = []
        # for l, w in self.w_lines:
        #     if segment.interval_crossed_by(l, self.left_x, self.right_x):
        #         up_w_lines.append((l, w))
        #         down_w_lines.append((l, w))
        #     elif l.above_interval(segment, self.left_x, self.right_x):
        #         up_w_lines.append((l, w))
        #     else:
        #         down_w_lines.append((l, w))
        #
        #
        # return Trapezoid(self.top_line, segment, self.left_x, self.right_x, up_w_lines, up), \
        #         Trapezoid(segment, self.bottom_line, self.left_x , self.right_x, down_w_lines, down)
        return Trapezoid(self.top_line, segment, self.left_x, self.right_x, self.w_lines[:], up), \
                Trapezoid(segment, self.bottom_line, self.left_x , self.right_x, self.w_lines[:], down)

    def vert_split(self, x):
        left = []
        right = []
        for p in self.points:
            if p[0] < x:
                left.append(p)
            else:
                right.append(p)
        return Trapezoid(self.top_line, self.bottom_line, self.left_x, x, self.w_lines, left), \
               Trapezoid(self.top_line, self.bottom_line, x, self.right_x, self.w_lines, right)


    def get_weight(self):
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
        elif (line.below_interval(self.top_line, self.left_x, self.right_x) and
                line.above_interval(self.bottom_line, self.left_x, self.right_x)):
                    return True
        elif line.interval_crossed_by_closed(self.top_line, self.left_x, self.right_x):
            return True
        else:
            return line.interval_crossed_by_closed(self.bottom_line, self.left_x, self.right_x)

    def get_vertices(self):
        def get_vertex(l, x):
            return (x, l.evaluate(x))
        vertices = [get_vertex(self.top_line, self.left_x), get_vertex(self.bottom_line, self.left_x),
                    get_vertex(self.top_line, self.right_x), get_vertex(self.bottom_line, self.right_x)]
        out_v = []
        for v in vertices:
            if abs(v[0]) != float("inf") and abs(v[1]) != float("inf"):
                out_v.append(v)
        return out_v


class Seidel_Tree:

    def is_active(self, node):
        return node.get_weight() > self.min_weight and len(node.get_points()) > self.min_p_count

    def register_leaf(self, node):
        if self.is_active(node):
            self.active_set.add(node)

    def unregister(self, node):
        if node in self.active_set:
            self.active_set.remove(node)

    def get_heaviest(self):
        return max(self.get_leaves(), key=lambda x: x.pt_count())

    def __init__(self, weighted_lines,  points=list(), min_weight=-1, min_p_count=-1):

        self.root = Trapezoid(w_lines=weighted_lines, points=points)
        self.lines = set()
        self.min_weight = min_weight
        self.min_p_count = min_p_count

        self.active_set = set()
        self.register_leaf(self.root)

    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def get_leaves(self) -> Trapezoid:
        stack = deque([self.root])

        all_traps = []
        while stack:
            curr_node = stack.pop()
            if curr_node.is_terminal():
                all_traps.append(curr_node)
            else:
                stack.append(curr_node.get_l_or_b())
                stack.append(curr_node.get_r_or_a())

        return all_traps

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

    def visualize_trapezoids(self, ax, min_x, max_x, min_y, max_y):
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
        all_traps = self.get_leaves()
        random.shuffle(all_traps)

        # ix = 0
        # curr_node = all_traps[6]
        # print(len(all_traps))
        # curr_node.visualize(ax, min_x, max_x, min_y, max_y, viz_lines=False, color=color_list[0])

        for curr_node in all_traps:
            curr_node.visualize(ax, min_x, max_x, min_y, max_y)
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
        self.unregister(terminal_node)
        self.register_leaf(left_term)
        self.register_leaf(right_term)

        if parent_node is None:
            self.root = Vertical_Node(x_sp, left_term, right_term)
        elif parent_node.get_l_or_b() == terminal_node:
            parent_node.set_l_or_b(Vertical_Node(x_sp, left_term, right_term))
        else:
            parent_node.set_r_or_a(Vertical_Node(x_sp, left_term, right_term))

        return left_term, right_term



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

        upper_zone = []
        lower_zone = []
        while node_stack:

            curr_node, parent_node, curr_segment = node_stack.pop()
            if curr_node.is_terminal():
                if self.is_active(curr_node):
                    upper, lower = curr_node.horz_split(curr_segment)
                    upper_zone.append((upper, parent_node))
                    lower_zone.append((lower, parent_node))
                    self.unregister(curr_node)
                    self.register_leaf(upper)
                    self.register_leaf(lower)

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

        def zone(self, line, curr_node=None):
            """
            Converts a line into many seperate segments.
            """

            curr_node = self.root if curr_node is None else curr_node
            full_segment = Segment(line, -float("inf"), float("inf"))
            node_stack = deque([(curr_node, None, full_segment)])

            while node_stack:

                curr_node, parent_node, curr_segment = node_stack.pop()
                if curr_node.is_terminal():
                    yield curr_node
                elif isinstance(curr_node, Segment_Node) and curr_node.segment.same_line(curr_segment):
                    continue
                elif curr_node.crosses(curr_segment):
                    if isinstance(curr_node, Segment_Node):
                        upper_split, lower_split, upper_right = curr_segment.split(curr_node.segment)
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

        def merge_traps(zone, key):
            line_matchings = {}
            for u_l, p in zone:
                if key(u_l) in line_matchings:
                    line_matchings[key(u_l)].append((u_l, p))
                else:
                    line_matchings[key(u_l)] = [(u_l, p)]

            for tp_lines in line_matchings:
                merge_traps = line_matchings[tp_lines][0][0]
                m_lines = set(merge_traps.w_lines)
                m_points = set(merge_traps.points)

                for (trap, parent_node) in line_matchings[tp_lines][1:]:
                    merge_traps.left_x = min(merge_traps.left_x, trap.left_x)
                    merge_traps.right_x = max(merge_traps.right_x, trap.right_x)
                    m_lines.update(trap.w_lines)
                    m_points.update(trap.points)
                    if parent_node.get_l_or_b() == trap:
                        parent_node.set_l_or_b(merge_traps)
                    else:
                        parent_node.set_r_or_a(merge_traps)
                merge_traps.w_lines = list(m_lines)
                merge_traps.points = list(m_points)

        if merge:
            merge_traps(upper_zone, lambda x: x.top_line)
            merge_traps(lower_zone, lambda x: x.bottom_line)



def to_line(p1, p2):
    return Line((p1[1] - p2[1]) / (p1[0] - p2[0]), p1[1] - (p1[1] - p2[1]) / (p1[0] - p2[0]) * p1[0])


def compute_cutting(test_set, weight_map, points, r) -> Seidel_Tree:
    total_weight = sum(weight_map[l] for l in test_set)
    min_weight = total_weight / r
    # Cutting size
    # min_pt_count =  len(points) / (12 * r * r) - 1
    lines = weighted_shuffle(test_set, [weight_map[l] / float(total_weight) for l in test_set])
    #random.shuffle(test_set)
    tree = Seidel_Tree([(l,weight_map[l]) for l in lines],
                       points,
                       min_weight=min_weight)
    for l in lines:
        if len(tree.active_set) == 0:
            return tree
        tree.insert_line(l)
    return tree

    #return tree


def plot_cutting_size(pts, l_s, h_s, k, r=4):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = []
        rand_pts = random.sample(pts, int(math.sqrt(i)))

        for p1, p2 in itertools.combinations(rand_pts, 2):
            test_set.append(to_line(p1, p2))
        random.shuffle(test_set)
        weight_map = {t: 1 for t in test_set}
        tree = compute_cutting(test_set, weight_map, pts, r)
        trap_count.append(len(tree.get_leaves()) / float(r * r))
        line_count.append(len(test_set))

    f, a = plt.subplots()
    a.scatter(line_count, trap_count)
    plt.show()

def plot_cutting_size3(pts, l_s, h_s, k, r=4):

    line_count = []
    trap_count = []
    for i in np.linspace(l_s, h_s, k):
        print(i)
        test_set = []
        for _ in range(int(i)):
            [p1, p2] = random.sample(pts, 2)
            test_set.append(to_line(p1, p2))
        weight_map = {t: 1 for t in test_set}
        tree = compute_cutting(test_set, weight_map, pts, r)
        trap_count.append(len(tree.get_leaves()) / float(r * r))
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
        trap_count.append(len(tree.get_leaves()) / float(r * r))
        line_count.append(len(test_set))

    f, a = plt.subplots()
    a.scatter(line_count, trap_count)
    plt.show()


# def find_bad_sequence():
#     global trigered
#
#     r = 4
#     while not trigered:
#         print("tested")
#
#         pts = [(random.random(), random.random()) for i in range(100)]
#         w_lines = []
#         rand_pts = random.sample(pts, 20)
#         lines = dual_cutting(pts, 2)
#         for l in lines:
#             w_lines.append((l, 1))
#         random.shuffle(w_lines)
#         total_weight = sum(w for l, w in w_lines)
#         min_weight = total_weight / r
#
#         tree = Seidel_Tree(w_lines, points=pts, min_weight=min_weight, min_p_count=-1)
#         for i in range(3):
#             new_line, _ = w_lines[i]
#             tree.insert_line(new_line)
#             if trigered:
#                 print(w_lines)
#                 print(pts)
#                 print(i)
#                 return







def remove_not_in_pts(pts, vertices):
    not_new = []
    for p1 in vertices:
        for p2 in pts:
            if approx_eq((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2, 0):
                not_new.append(p2)
                break
    return not_new



def visualize_cuttings(line_count, w_line_count, point_count, min_x, max_x, min_y, max_y, r, merge=False):
    pts = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
    w_lines = []
    rand_pts = random.sample(pts, int(math.sqrt(w_line_count)))
    for p1, p2 in itertools.combinations(rand_pts, 2):
        w_lines.append((to_line(p1, p2), 1))

    random.shuffle(w_lines)

    tree = Seidel_Tree(w_lines, points=pts, min_weight=-1, min_p_count=-1)
    matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
    for i in range(line_count):
        f, a = plt.subplots()
        a.set_ylim(min_y, max_y)
        a.set_xlim(min_x, max_x)
        new_line, _ = w_lines[i]

        tree.insert_line(new_line, merge=merge)
        tree.visualize_trapezoids(a, min_x, max_x, min_y, max_y)
        print(new_line)
        plt.show()


# def visualize_cuttings2(point_count, min_x, max_x, min_y, max_y, t, r2, merge=False):
#     pts = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
#     def one_gen():
#         while True:
#             yield 1
#
#     w_lines = list(zip(dual_cutting(pts, t), one_gen()))
#     print(len(w_lines))
#     print(w_lines)
#     random.shuffle(w_lines)
#
#     min_weight = len(w_lines) / r2 - 1
#     # Cutting size
#     min_pt_count = len(pts) / (18 * r2 * r2) - 1
#     tree = Seidel_Tree(w_lines, points=pts, min_weight=min_weight, min_p_count=-1)
#     #matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     for l in w_lines:
#         f, a = plt.subplots()
#         a.set_ylim(min_y, max_y)
#         a.set_xlim(min_x, max_x)
#         print("here")
#         tree.insert_line(l[0], merge=merge)
#         print("inserted")
#         tree.visualize_trapezoids(a, min_x, max_x, min_y, max_y)
#         plt.show()


# def visualize_cuttings3(point_count, min_x, max_x, min_y, max_y, t, r2, merge=False):
#     pts = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
#
#     lines = dual_cutting(pts, t)
#
#     random.shuffle(lines)
#
#     min_weight = len(lines) / r2 - 1
#     # Cutting size
#     min_pt_count = len(pts) / (18 * r2 * r2) - 1
#     tree = compute_cutting(lines, {l:1 for l in lines}, pts, r2)
#     matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
#     f, a = plt.subplots()
#     a.set_ylim(min_y, max_y)
#     a.set_xlim(min_x, max_x)
#
#     tree.visualize_trapezoids(a, min_x, max_x, min_y, max_y)
#     x, y = zip(*pts)
#
#
#     a.plot(list(x), list(y), 'x')
#     plt.show()


# pts = [(random.random(), random.random()) for k in range(10000)]
# plot_cutting_size(pts, 10, 10000, 10)

#find_bad_sequence()

#visualize_cuttings3(1000, -4, 4, -4, 4, 2, 16)

# point_count = 10000
# pts = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(point_count)]
# lines = []
# rand_pts = random.sample(pts, 20)
# for p1, p2 in itertools.combinations(rand_pts, 2):
#     lines.append(to_line(p1, p2))
# tree = compute_cutting(lines, {l:1 for l in lines}, pts, 16)


#error_plot(10000, 10, 400, 40, r=4, t=2)
# pts = [(4 * random.random() - 2, 4 * random.random() - 2) for i in range(10000)]
# lines = dual_cutting(pts, 2)
# print(len(lines))
# f, ax = plt.subplots()
# ax.set_ylim(-2, 2)
# ax.set_xlim(-2, 2)
# for l in lines:
#     ax.plot([-2, 2], [l.evaluate(-2), l.evaluate(2)])
# plt.show()


#

# plot_cutting_size(pts, 16, 4000, 80)
#plot_cutting_size2(pts, 16, 5000, 80)
# #compute_cutting()
#error_plot(pts, 10, 400, 40, r=2, z=50)


# tree = partition_tree(pts, 2, 20, z=1000)
# f, a = plt.subplots()
# # #
# tree.visualize(a, -5, 5, -5, 5, viz_lines=False, color = 'r')
# xs, ys = zip(*pts)
# a.scatter(xs, ys)
# plt.show()


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
#
# pts = [(random.random(), random.random()) for i in range(1000)]
# lines = test_set_lines(pts, 25)
# tree = compute_cutting(lines, {l: 1 for l in lines}, pts, 5)
# matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
# f, ax = plt.subplots()
# tree.visualize_trapezoids(ax, -1, 2, -1, 2)
# ax.set_xlim(0, 1)
# ax.set_ylim(0, 1)
# plt.show()
