from Cuttings import *
from collections import deque
import pydot

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

        e1 = Segment(self, self.xl, x_mid)
        e2 = Segment(self, x_mid, self.xr)
        return e1, e2

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
            return other.a == self.a and other.b == self.b and other.xl == self.xl and other.xr == self.xr


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

    def __init__(self, right, left):
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
    def __init__(self, x, right=None, left=None):
        super(Vertical_Node, self).__init__(right, left)
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
        super(Segment_Node, self).__init__(up, down)
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
                return self.segment.pt_eq_below(seg.right_vertex)
        else:
            if self.segment.contained(seg.left_vertex):
                return self.segment.below_closed_interval(seg, seg.xl, float("inf"))
            else:
                return self.segment.pt_eq_below(seg.left_vertex)

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

class Trapezoid(Node):
    def __repr__(self):
        return type(self).__name__ + "(**" + pprint.pformat(vars(self), indent=4, width=1) + ")"

    def __init__(self, top_line=None, bottom_line=None, left_x=-float("inf"), right_x=float("inf"), lines=list(), weights=list()):
        self.top_line = Line(0, float("inf")) if top_line is None else top_line

        self.bottom_line = Line(0, float("-inf")) if bottom_line is None else bottom_line
        self.left_x = left_x
        self.right_x = right_x
        self.lines = lines
        self.weights = weights
        self.weight = sum(weights)

    def inside(self, pt):
        """
        Determines if the point is inside
        :param pt:
        :return:
        """
        pass

    def is_terminal(self):
        return True

    def horz_split(self, segment):
        """
        This can only cut in two, because we always do vertical splits first.
        :param segment:
        :return:
        """
        bottom = []
        top = []
        tw = []
        bw = []
        for l, w in zip(self.lines, self.weights):
            #if the line crosses the segment inside of the trapezoid
            if segment.interval_crossed_by(l, self.left_x, self.right_x):
                bottom.append(l)
                bw.append(w)
                top.append(l)
                tw.append(w)
            #if the line is above the segment inside the trapezoid
            elif l.above_closed_interval(segment, segment.xl, segment.xr):
                top.append(l)
                tw.append(w)
            else:
                bottom.append(l)
                bw.append(w)
        return Trapezoid(self.top_line, segment, self.left_x, self.right_x, top, tw), \
                Trapezoid(segment, self.bottom_line, self.left_x , self.right_x, bottom, bw)


    def vert_split(self, x):
        left = []
        right = []
        lw = []
        rw = []
        ty = self.top_line.evaluate(x)
        by = self.bottom_line.evaluate(x)
        for l, w in zip(self.lines, self.weights):
            ly = l.evaluate(x)
            if approx_above(by, ly) and approx_above(ly, ty):
                left.append(l)
                right.append(l)
                lw.append(w)
                rw.append(w)
            else:
                wedge_point = self.top_line.x_intercept(self.bottom_line)
                intercepts = [l.x_intercept(self.top_line), l.x_intercept(self.bottom_line)]
                if wedge_point < x:
                    intercepts = [x for x in intercepts if x > wedge_point]
                else:
                    intercepts = [x for x in intercepts if x < wedge_point]
                if intercepts[0] < x:
                    left.append(l)
                    lw.append(w)
                else:
                    right.append(l)
                    rw.append(w)
        return Trapezoid(self.top_line, self.bottom_line, self.left_x, x, left, lw), \
                Trapezoid(self.top_line, self.bottom_line, x, self.right_x, right, rw)

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


class Seidel_Tree:
    def  __init__(self, crossing_lines=list(), weights=None):
        if weights is None:
            weights = [1 for i in range(len(crossing_lines))]
        self.root = Trapezoid(lines=crossing_lines, weights=weights)

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

        def gen_label(n):
            if isinstance(n, Segment_Node):
                name = "Line(%f, %f, %f, %f)" % (n.segment.a, n.segment.b, n.segment.xl, n.segment.xr)
            elif isinstance(n, Trapezoid):
                name = "T(%d)"%(n.weight,)
            elif isinstance(n, Vertical_Node):
                name = "V(%f)" % (n.x,)
            return name
        stack = deque([(self.root, 0)])
        graph = pydot.Dot(graph_type='graph')
        while stack:
            curr_node, nid = stack.pop()
            node = pydot.Node("%s_%d"% (toOutputString(curr_node), nid), label=gen_label(curr_node))
            graph.add_node(node)
            if not curr_node.is_terminal():
                if curr_node.get_r_or_a() is not None:
                    edge = pydot.Edge("%s_%d"% (toOutputString(curr_node), nid), "%s_%d" %
                                 (toOutputString(curr_node.get_r_or_a()), 2 * nid + 1))
                    graph.add_edge(edge)
                    stack.append((curr_node.get_r_or_a(), 2 * nid + 1))
                if curr_node.get_l_or_b() is not None:
                    edge = pydot.Edge("%s_%d" % (toOutputString(curr_node), nid), "%s_%d" %
                                  (toOutputString(curr_node.get_l_or_b()), 2 * nid + 2))

                    graph.add_edge(edge)
                    stack.append((curr_node.get_l_or_b(), 2 * nid + 2))
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
                    print(segments)
                    top_y = max_val(segments, curr_node.x)
                    bottom_y = min_val(segments, curr_node.x)
                    if min_x < curr_node.x < max_x:
                        ax.vlines(curr_node.x, max(bottom_y, min_y), min(top_y, max_y))
                if curr_node.get_r_or_a() is not None:
                    stack.append((curr_node.get_r_or_a(), segments))
                if curr_node.get_l_or_b() is not None:
                    stack.append((curr_node.get_l_or_b(), segments))

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

    def left_child(self, segment, curr_node):
        if curr_node.is_terminal():
            return None
        elif curr_node.get_l_or_b() and curr_node.get_l_or_b().crosses(segment):
            return curr_node.get_l_or_b()
        else:
            return None

    def right_child(self, segment, curr_node):
        if curr_node.is_terminal():
            return None
        elif curr_node.get_r_or_a() and curr_node.get_r_or_a().crosses(segment):
            return curr_node.get_r_or_a()
        else:
            return None

    def line_to_segments(self, line, curr_node=None):
        """
        Converts a line into many seperate segments.
        """
        curr_node = self.root if curr_node is None else curr_node
        full_segment = Segment(line, -float("inf"), float("inf"))
        node_stack = deque()
        x_vals = [-float("inf")]
        while True:
            if curr_node is not None:
                node_stack.append(curr_node)
                curr_node = self.left_child(full_segment, curr_node)
            else:
                if not node_stack:
                    break
                curr_node = node_stack.pop()
                if isinstance(curr_node, Segment_Node): #and curr_node.segment.crossed_by(line):
                    x_vals.append(curr_node.segment.x_intercept(full_segment))
                curr_node = self.right_child(full_segment, curr_node)

        # x_vals should be in increasing order
        x_vals.append(float("inf"))
        segments = []
        print(x_vals)
        x_vals = sorted(x_vals)
        for xl, xr in zip(x_vals[:-1], x_vals[1:]):
            segments.append(Segment(line, xl, xr))
        return segments

    def insert_splitter(self, line, right, curr_node=None, min_weight=-1):
        """
        :return:
        """
        curr_node = self.root if curr_node is None else curr_node
        terminal_node, parent_node = self.locate_segment_end(line, right, curr_node=curr_node)

        if right:
            x_sp = line.xr
        else:
            x_sp = line.xl
        if isinstance(parent_node, Vertical_Node) and approx_eq(parent_node.x, x_sp):
            return

        if terminal_node.getWeight() <= min_weight:
            return
        left_term, right_term = terminal_node.vert_split(x_sp)
        if parent_node is None:
            self.root = Vertical_Node(x_sp, left_term, right_term)
        elif parent_node.get_l_or_b() == terminal_node:
            parent_node.set_l_or_b(Vertical_Node(x_sp, left_term, right_term))
        else:
            parent_node.set_r_or_a(Vertical_Node(x_sp, left_term, right_term))

    def insert_segment(self, seg, curr_node=None, min_weight=-1):
        """
        The segment is assumed to not cross other line segments.
        :return:
        """
        curr_node = self.root if curr_node is None else curr_node
        for (n, p) in self.locate_segment(seg, curr_node=curr_node):
            if n.getWeight() <= min_weight:
                continue
            upper, lower = n.horz_split(seg)
            new_s = Segment_Node(seg, upper, lower)
            if p is None:
                self.root = new_s
            elif p.get_l_or_b() == n:
                p.set_l_or_b(new_s)
            else:
                p.set_r_or_a(new_s)

    def insert_full_segment(self, seg, curr_node=None, min_weight=-1):
        """
        The segment is assumed to not cross other line segments.
        :return:
        """
        curr_node = self.root if curr_node is None else curr_node
        if seg.xl != float("-inf"):
            self.insert_splitter(seg, False, curr_node, min_weight=min_weight)
        if seg.xr != float("inf"):
            self.insert_splitter(seg, True, curr_node, min_weight=min_weight)
        self.insert_segment(seg, curr_node, min_weight=min_weight)

    def insert_line(self, line, curr_node=None, min_weight=-1):
        curr_node = self.root if curr_node is None else curr_node
        segments = self.line_to_segments(line, curr_node=curr_node)
        for seg in segments:
            if approx_eq(seg.xl, seg.xr):
                continue
            print(seg)

            self.insert_full_segment(seg, curr_node=curr_node, min_weight=min_weight)


if __name__ == "__main__":
    #arrange = Arrangement([random_box_line() for i in range(10)])
    #assume the domain is a 1 by 1 square.
    #construct the lines so they extend from the left side to the right side.
    max_x = 2
    min_x = -2
    tree = Seidel_Tree([random_box_line() for i in range(100)], [1 for i in range(100)])
    matplotlib.rcParams['figure.figsize'] = [20.0, 20.0]
    for i in range(6):
        #f, (a1, a2) = plt.subplots(2, 1)
        new_line = random_box_line()
        #print(new_line)
        #arrange.insert(new_line)
        tree.insert_line(new_line, min_weight=10)

        #print(tree)
    tree.visualize("testing")
    f, a = plt.subplots()
    tree.visualize_arrangement(a, -5, 5, -5, 5)
    plt.show()