import random
import itertools
import SeidelTree as St
import Cuttings as Ct
import math

def to_dual_line(pt):
    return St.Line(-pt[0], pt[1])

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
    dual_lines = []
    for p in pts:
        dual_lines.append(to_dual_line(p))
    random.shuffle(dual_lines)
    tree = St.compute_cutting(dual_lines, {l:1 for l in dual_lines}, [], r)
    vertices = []
    for trap in tree.get_leaves():
        vertices.extend(trap.get_vertices())
    test_set = []
    vertices = deduplicate_points(vertices)

    test_set = []
    for v in vertices:
        test_set.append(to_dual_line(v))
    return test_set


def test_set_dual(pts, t):
    t_tmp = t
    while True:
        lines = dual_cutting(pts, max(int(math.sqrt(t_tmp / 8) + 1), 2))
        if len(lines) < t:
            t_tmp = 2 * t_tmp
            print("doubling")
        else:
            break
    random.shuffle(lines)
    return lines[:t]


def test_set_points(pts, t):
    test_set = []
    rand_pts = random.sample(pts, int(math.sqrt(2 * t) + 1))
    for p1, p2 in itertools.combinations(rand_pts, 2):
        if len(test_set) >= t:
            break
        test_set.append(St.to_line(p1, p2))
    return test_set


def test_set_lines(pts, t):
    test_set = []
    rand_pts = random.sample(pts, 2 * t)

    for p1, p2 in zip(rand_pts[:t], rand_pts[t:]):
        test_set.append(St.to_line(p1, p2))
    return test_set

