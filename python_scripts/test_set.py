import random
import itertools
import SeidelTree as St
import PolyTree as Pt
import Cuttings as Ct
import math





def deduplicate_points(pts):
    if not pts:
        return []
    s_pts_x = sorted(pts, key=lambda x: x[0])

    curr_p = s_pts_x[0]
    groups = [[curr_p]]
    for p in s_pts_x[1:]:
        if Ct.approx_eq(curr_p[0], p[0]):
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
            if not Ct.approx_eq(c_p[1], p[1]):
                not_duplicate.append(p)
                c_p = p

    return not_duplicate


def dual_cutting(pts, r):
    #print(r)
    dual_lines = []
    for p in pts:
        dual_lines.append(Ct.to_dual_line(p))
    #print("dual_lines=%d ,r = %d"%(len(dual_lines), r))
    tree = Pt.compute_cutting(dual_lines, {l:1 for l in dual_lines}, [], r)
    #print("computed")
    vertices = []
    #print(len(tree.get_leaves()))
    for trap in tree.get_leaves():
        vertices.extend(trap.get_vertices())
    print(len(vertices))
    vertices = deduplicate_points(vertices)
    #print("got here")
    test_set = []
    for v in vertices:
        test_set.append(to_dual_line(v))
    return test_set, tree


def test_set_dual_exact_t(pts, t, c=1):
    internal_pts = random.sample(pts, int(.5 + min(c * math.sqrt(t) * math.log(t), len(pts))))
    return dual_cutting(internal_pts, max(int((math.sqrt(t)) + .1), 1))[0]


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