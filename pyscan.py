from libpyscan import *
import random
import bisect
import itertools
import collections

def to_weighted(points):
    return [WPoint(1.0, pt[0], pt[1], 1.0) for pt in points]

def trajectories_to_flux(trajectories):
    return [t[0] for t in trajectories], [t[-1] for t in trajectories]

def trajectories_to_labels(trajectories):
    return itertools.chain.from_iterable([[LPoint(label, 1.0, pt[0], pt[1], 1.0) for pt in traj]
            for label, traj in zip(range(len(trajectories)), trajectories)])

def my_sample(samp, count):
    return random.sample(samp, min(len(samp), int(count + .5)))

def evaluate_range(range, mp, bp, disc_f):

    if not mp and not bp:
        return evaluate(disc_f, 0, 0, 0, 0)
    elif not mp:
        pt_obj = bp[0]
    else:
        pt_obj = mp[0]

    if isinstance(range, Disk):
        if isinstance(pt_obj, LPoint):
            return evaluate_disk_labeled(range, mp, bp, disc_f)
        elif isinstance(pt_obj, WPoint) or isinstance(pt_obj, tuple):
            return evaluate_disk(range, mp, bp, disc_f)
    elif isinstance(range, Halfplane):
        if isinstance(pt_obj, LPoint):
            return evaluate_halfplane_labeled(range, mp, bp, disc_f)
        elif isinstance(pt_obj, WPoint) or isinstance(pt_obj, tuple):
            return evaluate_halfplane(range, mp, bp, disc_f)
    elif isinstance(range, Rectangle):
        if isinstance(pt_obj, LPoint):
            return evaluate_rectangle_labeled(range, mp, bp, disc_f)
        elif isinstance(pt_obj, WPoint) or isinstance(pt_obj, tuple):
            return evaluate_rectangle(range, mp, bp, disc_f)
    raise ValueError()

def evaluate_range_trajectory(range, mp, bp, disc_f):
    """
    Evaluates this range to compute the total discrepancy.
    :param range:
    :param mp:
    :param bp:
    :param disc_f:
    :return:
    """
    if not mp and not bp:
        return evaluate(disc_f, 0, 0, 0, 0)
    if isinstance(range, Disk):
        return evaluate_disk_trajectory(range, mp, bp, disc_f)
    elif isinstance(range, Halfplane):
        return evaluate_halfplane_trajectory(range, mp, bp, disc_f)
    elif isinstance(range, Rectangle):
        return evaluate_rectangle_trajectory(range, mp, bp, disc_f)
    else:
        raise ValueError()


def split_set(pts, rate):
    """
    Divides the point set into two random sets where one contains approximately rate * len(pts) number of
    points and the other has (1 - rate) * len(pts)
    :param pts:
    :param rate:
    :return:
    """
    red_set = my_sample(pts, len(pts) * rate)
    red_set_set = set(red_set)
    return red_set, [item for item in pts if item not in red_set_set]




"""
This takes a scanning function and two point sets and then computes a planted region that contains some fraction r
of the points with some tolerance.
This then computes a red and blue set of points based on the planted region.

Inside the region q fraction of points are red.
Outside the region p fraction of points are red
"""
def plant_region(points, r, p, q, eps, scan_f):

    while True:

        net_size = int(1 / min(r, eps) + 1)
        s_size = int(1 / min(r, eps)** 2 + 1)

        net = my_sample(points, net_size)
        sample = my_sample(points, s_size)
        disc = size_region(r)
        reg, _ = scan_f(net, sample, [], disc)

        val = evaluate_range(reg, points, [], disc)
        if 1 - val < eps:
            break

    in_region = []
    out_region = []
    for pt in pts:
        if reg.contains(pt):
            in_region.append(pt)
        else:
            out_region.append(pt)

    red_in, blue_in = split_set(in_region, q)

    red_out, blue_out = split_set(out_region, p)
    return red_in + red_out, blue_in + blue_out, reg


def plant_full_square(trajectories, r, p, q, disc, max_count=32):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    if not trajectories:
        return None
    count = 0
    traj = []
    while not traj:
        count+= 1
        traj = random.choice(trajectories)
        if count > 10 :
            return None

    seed_pt = random.choice(traj)
    upper_bound = 1.0
    lower_bound = 0.0
    num = 0
    while num < max_count:
        size = (upper_bound + lower_bound) / 2
        reg = Rectangle(seed_pt[0] + size / 2, seed_pt[1] + size / 2, seed_pt[0] - size / 2, seed_pt[1] - size / 2)
        count = sum(1 for traj in trajectories if reg.intersects_trajectory(traj))
        if abs(count - r * len(trajectories)) <= 2:
            break
        if count - r * len(trajectories) > 0:
            upper_bound = size
        else:
            lower_bound = size
        num += 1


    inside_rect = [traj for traj in trajectories if reg.intersects_trajectory(traj)]
    red_in, blue_in = split_set([tuple(traj) for traj in inside_rect], q)
    outside_rect = [traj for traj in trajectories if not reg.intersects_trajectory(traj)]
    red_out, blue_out = split_set([tuple(traj) for traj in outside_rect], p)

    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, reg, diff




def plant_full_halfplane(trajectories, r, p, q, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    def min_distance(pts, direc):
        return min([direc[0] * pt[0] + direc[1] * pt[1] for pt in pts])
    rand_direc = (random.gauss(0,1), random.gauss(0, 1))
    trajectories = sorted(trajectories, key=lambda x: min_distance(x, rand_direc))

    pt = min(trajectories[int(r * len(trajectories))], key=lambda pt: rand_direc[0] * pt[0] + rand_direc[1] * pt[1])
    lc = pt[0] * rand_direc[0] + pt[1] * rand_direc[1]
    plant_region = Halfplane(Point(rand_direc[0], rand_direc[1], -lc))

    inside_plane = trajectories[:int(r * len(trajectories))]
    outside_plane = trajectories[int(r * len(trajectories)):]
    red_in, blue_in = split_set([tuple(traj) for traj in inside_plane], q)
    red_out, blue_out = split_set([tuple(traj) for traj in outside_plane], p)

    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, plant_region, diff


def plant_partial_halfplane(trajectories, r, p, q, eps, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    def min_distance(pt, direc):
        return direc[0] * pt[0] + direc[1] * pt[1]
    #so it always points down
    rand_direc = (random.gauss(0,1), -abs(random.gauss(0, 1)))

    trajectory_obj = [Trajectory(pts) for pts in trajectories]
    all_pts = uniform_sample(trajectory_obj, int(1 / eps ** 2 + 1), False)
    all_pts = sorted(all_pts, key=lambda x: min_distance(x, rand_direc))

    pt = all_pts[int((1 - r) * len(all_pts))]
    lc = pt[0] * rand_direc[0] + pt[1] * rand_direc[1]
    plant_region = Halfplane(Point(rand_direc[0], rand_direc[1], -lc))

    inside_disk = [traj for traj in trajectory_obj if plant_region.intersects_trajectory(traj)]
    outside_disk = [traj for traj in trajectory_obj if not plant_region.intersects_trajectory(traj)]

    red_in, blue_in = split_set(inside_disk, q)
    red_out, blue_out = split_set(outside_disk, p)
    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, plant_region, diff


def plant_full_disk(trajectories, r, p, q, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    origin = random.choice(list(itertools.chain.from_iterable(trajectories)))
    trajectory_obj = [Trajectory(pts) for pts in trajectories]

    trajectories = sorted(trajectory_obj, key=lambda el: el.point_dist(origin))

    inside_disk = trajectories[:int(r * len(trajectories))]
    outside_disk = trajectories[int(r * len(trajectories)):]

    max_disk = Disk(origin[0], origin[1], inside_disk[-1].point_dist(origin))

    red_in, blue_in = split_set([tuple(traj.get_pts()) for traj in inside_disk], q)
    red_out, blue_out = split_set([tuple(traj.get_pts()) for traj in outside_disk], p)
    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, max_disk, diff


def plant_partial_disk(trajectories, r, p, q, eps, disc):
    """
    Choose a point at random from a trajectory and then expand outward from there.
    :param trajectories this consists of lists of lists of points where the points are type Pyscan.Point
    :param r:
    :param p:
    :param q:
    :return:
    """
    #Compute fraction of mass inside disk
    #For each segment inside disk add mass.
    # For each segment partly inside disk compute overlap and add to mass.
    #We can do bisection on this amount since it is a monotonic function to compute within eps fraction.

    origin_tuple = random.choice(list(itertools.chain.from_iterable(trajectories)))
    origin = Point(origin_tuple[0], origin_tuple[1], 1.0)
    trajectory_obj = [Trajectory(pts) for pts in trajectories]
    all_pts = uniform_sample(trajectory_obj, int(1 / eps ** 2 + 1), False)
    trajectories = sorted(all_pts, key=lambda x: x.dist(origin))

    disk_boundary = trajectories[int(r * len(all_pts))]
    max_disk = Disk(origin[0], origin[1], origin.dist(disk_boundary))
    inside_disk = [traj for traj in trajectory_obj if max_disk.intersects_trajectory(traj)]
    outside_disk = [traj for traj in trajectory_obj if not max_disk.intersects_trajectory(traj)]

    red_in, blue_in = split_set(inside_disk, q)
    red_out, blue_out = split_set(outside_disk, p)
    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, max_disk, diff

"""
Create a distribution of the null distribution to measure the significance of the region.
"""
def distribution(points, p, scan_f, n, s, disc=DISC):
    while True:
        red, blue = split_set(points, p)

        net_set = my_sample(red, min(len(red), n)) + my_sample(blue, min(len(blue), n))
        m_sample = my_sample(red, min(len(red), s), red)
        b_sample = my_sample(blue, min(len(blue), s), blue)
        reg, val = scan_f(net_set, m_sample, b_sample, disc)

        yield val


def null_cdf(observations):
    """
    Generates a cdf object using a certain number of observations.
    :param observations
    :return:
    """
    values = sorted(observations)
    prob = [x / float(len(observations)) for x in range(len(observations))]
    return values, prob


def random_rect(points, r):
    def bBox(*bb):
        return (min(*map(lambda pt: pt[0], bb)),
                max(*map(lambda pt: pt[0], bb)),
                min(*map(lambda pt: pt[1], bb)),
                max(*map(lambda pt: pt[1], bb)))

    while True:
        [pt1, pt2] = random.sample(points, 2)

        if random.randint(0, 1) == 0:
            (minX, maxX, minY, maxY) = bBox(pt1, pt2)
            fx = lambda pt: pt[0]
            fy = lambda pt: pt[0]
        else:
            fx = lambda pt: pt[0]
            fy = lambda pt: pt[0]
            (minY, maxY, minX, maxX) = bBox(pt1, pt2)

        if random.randint(0, 1) == 0:
            sf = lambda pt: fx(pt)
        else:
            sf = lambda pt: -fx(pt)
            minX, maxX = maxX, minX

        subPoints = [pt for pt in points if minY <= fy(pt) and fy(pt) <= maxY]
        # Make a pass through the data and find the lower and upper bound.
        inB = [pt for pt in subPoints if minX <= fx(pt) and fx(pt) <= maxX]
        if random.randint(0, 1):  # random.randint(0, 1) == 0:
            ab = [pt for pt in subPoints if maxX < fx(pt)]
            ab.sort(key=lambda p: -sf(p))
        else:
            ab = [pt for pt in subPoints if fx(pt) < minX]
            ab.sort(key=lambda p: sf(p))
        while len(ab) > 0:
            if len(inB) > int(r * len(points) + .5):
                break
            if len(inB) == int(r * len(points) + .5):
                (lx, ux, ly, uy) = bBox(*inB)
                return Rectangle(ux, uy, lx, ly, 0.0)
            el = ab.pop()
            inB.append(el)
            while len(ab) > 0 and sf(el) == sf(ab[-1]):
                el = ab.pop()
                inB.append(el)


def plant_rectangle(pts, r, p, q):
    """
    Create a set of red and blue points with a random rectangle planted containing r fraction of the points.
    :param p:
    :param q:
    :return:
    """
    rect = random_rect(pts, r)

    inside_rect = []
    outside_rect = []
    for pt in pts:
        if rect.contains(pt):
            inside_rect.append(pt)
        else:
            outside_rect.append(pt)

    red_in = split_set(inside_rect, q)
    blue_in = split_set(inside_rect, q)
    red_out = split_set(outside_rect, p)
    blue_out = split_set(outside_rect, p)
    return red_in + red_out, blue_in + blue_out, rect


def plant_disk(pts, r, p, q):
    """
    Create a set of red and blue points with a random disk planted containing r fraction of the points.
    :param p:
    :param q:
    :return:
    """

    selected = my_sample(pts, 1)
    if not selected:
        return None
    origin = selected[0]

    ordered = sorted(pts, key=lambda x: x.dist(origin))
    disk_boundary = ordered[int(r * len(ordered))]
    max_disk = Disk(origin[0], origin[1], origin.dist(disk_boundary))

    inside_disk = ordered[:int(r * len(ordered))]
    outside_disk = ordered[int(r * len(ordered)):]

    red_in = split_set(inside_disk, q)
    blue_in = split_set(inside_disk, q)
    red_out = split_set(outside_disk, p)
    blue_out = split_set(outside_disk, p)
    return red_in + red_out, blue_in + blue_out, max_disk


def plant_halfplane(pts, r, p, q):
    """
    Create a set of red and blue points with a random halfplane planted containing r fraction of the points.
    :param p:
    :param q:
    :return:
    """

    def min_distance(pt, direc):
        return direc[0] * pt[0] + direc[1] * pt[1]
    rand_direc = (random.gauss(0, 1), abs(random.gauss(0, 1)))

    ordered = sorted(pts, key=lambda x: min_distance(x, rand_direc))
    inside_halfplane = ordered[:int(r * len(ordered))]
    outside_halfplane = ordered[int(r * len(ordered)):]

    pt = ordered[int((1 - r) * len(ordered))]
    lc = pt[0] * rand_direc[0] + pt[1] * rand_direc[1]
    plant_region = Halfplane(Point(rand_direc[0], rand_direc[1], -lc))

    red_in = split_set(inside_halfplane, q)
    blue_in = split_set(inside_halfplane, q)
    red_out = split_set(outside_halfplane, p)
    blue_out = split_set(outside_halfplane, p)
    return red_in + red_out, blue_in + blue_out, plant_region


def plant_partial_rectangle(trajectories, r, p, q, eps, disc):

    trajectory_obj = [Trajectory(pts) for pts in trajectories]
    all_pts = uniform_sample(trajectory_obj, int(1 / eps ** 2 + 1), False)

    _, _, rect = plant_rectangle(all_pts, r, p, q)
    inside_rect = [traj for traj in trajectory_obj if rect.intersects_trajectory(traj)]
    outside_rect = [traj for traj in trajectory_obj if not rect.intersects_trajectory(traj)]
    red_in, blue_in = split_set(inside_rect, q)
    red_out, blue_out = split_set(outside_rect, p)

    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, rect, diff



"""
This plants a region where every trajectory:
Completely outside or inside of the region has an endpoint chosen at random.
Every trajectory with one endpoint inside the region has an endpoint chosen inside
with probability q (exactly q fraction have one endpoint in the region)

r controls how many points the region contains.

"""
def paired_plant_region(traj_start, traj_end, r, q, region_plant_f):

    _, _, reg = region_plant_f(traj_start + traj_end, r, .5, q)

    flux_region = []
    out_region = []
    for st_pt, end_pt in zip(traj_start, traj_end):
        if reg.contains(st_pt) and not reg.contains(end_pt):
            flux_region.append((st_pt, end_pt))
        elif reg.contains(end_pt) and not reg.contains(st_pt):
            flux_region.append((end_pt, st_pt))
        else:
            out_region.append((st_pt, end_pt))

    q_fraction, remainder = split_set(flux_region, q)
    remainder = [(ep, sp) for (sp, ep) in remainder]

    q_fraction_o, remainder_o = split_set(out_region, .5)
    remainder = [(ep, sp) for (sp, ep) in remainder]
    red, blue = zip(*(q_fraction + q_fraction_o + remainder + remainder_o))
    return red, blue, reg

