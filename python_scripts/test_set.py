import random
import itertools
import seidel_tree as St
import poly_tree as Pt
import geometric as Ct
import math

def dual_cutting(pts, r):
    #print(r)
    dual_lines = []
    for p in pts:
        dual_lines.append(Ct.to_dual_line(p))
    #print("dual_lines=%d ,r = %d"%(len(dual_lines), r))
    tree = Pt.PolyTree2(points=pts, lines=dual_lines)
    tree.cutting_r(r, {l: 1 for l in dual_lines})
    #print("computed")
    vertices = []
    #print(len(tree.get_leaves()))
    for trap in tree.get_leaves():
        vertices.extend(trap.get_vertices())
    #print(len(vertices))
    vertices = Ct.deduplicate_points(vertices)
    #print("got here")
    test_set = []
    for v in vertices:
        test_set.append(Ct.to_dual_line(v))
    return test_set, tree


def test_set_dual_exact_t(pts, t, c=1):
    internal_pts = random.sample(pts, int(.5 + min(c * math.sqrt(t) * math.log(t), len(pts))))
    return dual_cutting(internal_pts, max(int((math.sqrt(t)) + .1), 1))[0]


def maintain_test_set(tree, t):
    class smart_dict(dict):
        def __missing__(self, key):
            return 1
    tree.cutting_r(max(int((math.sqrt(t)) + .1), 2), smart_dict())
    vertices = []
    #print(len(tree.get_leaves()))
    for trap in tree.get_leaves():
        vertices.extend(trap.get_vertices())
    #print(len(vertices))
    vertices = Ct.deduplicate_points(vertices)
    #print("got here")
    test_set = []
    for v in vertices:
        test_set.append(Ct.to_dual_line(v))
    return test_set


def test_dual_tree(pts, t, c=1):
    internal_pts = random.sample(pts, int(.5 + min(c * math.sqrt(t) * math.log(t), len(pts))))
    return dual_cutting(internal_pts, max(int((math.sqrt(t)) + .1), 1))



def test_set_dual(pts, t, c = 1):
    internal_pts = random.sample(pts, int(.5 + min(c * math.sqrt(t) * math.log(t), len(pts))))
    t_tmp = t
    while True:
        lines = dual_cutting(internal_pts, max(int(math.sqrt(t_tmp) + 1), 2))[0]
        if len(lines) < t:
            t_tmp = 2 * t_tmp
            print("doubling")
        else:
            break
    random.shuffle(lines)
    return lines[:t]


def test_set_points(pts, t):
    test_set = []
    rand_pts = random.sample(pts, min(int(math.sqrt(2 * t) + 1), len(pts)))
    for p1, p2 in itertools.combinations(rand_pts, 2):
        if len(test_set) >= t:
            break
        test_set.append(St.to_line(p1, p2))
    return test_set


def test_set_lines(pts, t):
    test_set = []
    t = min(t, int(len(pts)/2))
    rand_pts = random.sample(pts, 2 *t)

    for p1, p2 in zip(rand_pts[:t], rand_pts[t:]):
        test_set.append(St.to_line(p1, p2))
    return test_set

#
# import matplotlib.pyplot as plt
#
# pts = [(random.random(), random.random()) for i in range(10000)]
#
# lines = test_set_dual_exact_t(pts, 100)
# print(len(lines))
# f, ax = plt.subplots()
# for l in lines:
#     l.visualize(ax, -1, 1)
# ax.set_xlim([-1, 1])
# ax.set_ylim([-1, 1])
# plt.show()