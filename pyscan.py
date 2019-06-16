from libpyscan import *
import libpyscan as lp
import random
import itertools
import math



def bounding_box(pts):
    if not pts:
        return None
    min_x = pts[0][0]
    max_x = pts[0][0]
    min_y = pts[0][1]
    max_y = pts[0][1]
    for pt in pts:
        min_x = min(min_x, pt[0])
        max_x = max(max_x, pt[0])
        min_y = min(min_y, pt[1])
        max_y = max(max_y, pt[1])
    return ((min_x, min_y), (max_x, max_y))


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
    """
    Evaluates this range to compute the total discrepancy over a set of points.

    :param range: Some arbitrary range.
    :param mp: measured set
    :param bp: baseline set
    :param disc_f: Discrepancy function.
    :return: Discrepancy function value.
    """
    if not mp and not bp:
        return evaluate(disc_f, 0, 0, 0, 0)
    elif not mp:
        pt_obj = bp[0]
    else:
        pt_obj = mp[0]

    if isinstance(pt_obj, LPoint):
        if isinstance(range, Disk):
            return evaluate_disk_labeled(range, mp, bp, disc_f)
        elif isinstance(range, Halfplane):
            return evaluate_disk_labeled(range, mp, bp, disc_f)
        elif isinstance(range, Rectangle):
            return evaluate_rectangle_labeled(range, mp, bp, disc_f)
    elif isinstance(pt_obj, WPoint):
        if isinstance(range, Disk):
            return evaluate_disk(range, mp, bp, disc_f)
        elif isinstance(range, Halfplane):
            return evaluate_halfplane(range, mp, bp, disc_f)
        elif isinstance(range, Rectangle):
            return evaluate_rectangle(range, mp, bp, disc_f)
    raise ValueError()


def evaluate_range_trajectory(range, mp, bp, disc_f):
    """
    Evaluates this range to compute the total discrepancy over a set of trajectories.

    :param range: Some arbitrary range.
    :param mp: measured set
    :param bp: baseline set
    :param disc_f: Discrepancy function.
    :return: Discrepancy function value.
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
    red_set = my_sample(pts, len(pts) * rate)
    red_set_set = set(red_set)
    return red_set, [item for item in pts if item not in red_set_set]


def plant_region(points, r, p, q, eps, scan_f):
    """
    This takes a scanning function and two point sets and then computes a planted region that contains some fraction r
    of the points with some tolerance.
    This then computes a red and blue set of points based on the planted region.

    Inside the region q fraction of points are red.
    Outside the region p fraction of points are red

    :param points: List of points.
    :param r: Fraction of points inside the region.
    :param p: Fraction of points outside the region that are red.
    :param q: Fraction of points inside the region that are red.
    :param eps: Difference between r and the fraction of points inside the region.
    :param scan_f: The scan function to use (ex max_disk)
    :return: red set, blue set, region planted.
    """
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


def plant_full_square(trajectories, r, p, q, max_count=32):
    """
    Choose a point at random from a trajectory and then expand a square out from this point till this region contains
    r fraction of all the trajectories.

    :param trajectories: List of list of points.
    :param r: double between 0 and 1
    :param p: double between 0 and 1
    :param q: double between 0 and 1
    :param disc: The discrepancy function to evaluate exactly on this region.
    :param max_count: The maximum number of times we will attempt to find the right sized region.
    :return: red set of trajectories, blue set of trajectories, the planted region, and the exact partial discrepancy.
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
        count = sum(1 for traj in trajectories if reg.intersects_trajectory(Trajectory(traj)))
        if abs(count - r * len(trajectories)) <= 2:
            break
        if count - r * len(trajectories) > 0:
            upper_bound = size
        else:
            lower_bound = size
        num += 1


    inside_rect = [traj for traj in trajectories if reg.intersects_trajectory(Trajectory(traj))]
    red_in, blue_in = split_set([tuple(traj) for traj in inside_rect], q)
    outside_rect = [traj for traj in trajectories if not reg.intersects_trajectory(Trajectory(traj))]
    red_out, blue_out = split_set([tuple(traj) for traj in outside_rect], p)

    diff = evaluate(disc, len(red_in), len(red_in) + len(red_out), len(blue_in), len(blue_in) + len(blue_out))
    return red_in + red_out, blue_in + blue_out, reg



def plant_full_halfplane(trajectories, r, p, q):
    """
    Choose a random direction and then finds a halfplane with this normal containing r fraction of the total trajectories.

    :param trajectories: List of list of points.
    :param r: Fraction of trajectories in region.
    :param p: Fraction of red trajectories outside of region.
    :param q: Fraction of red trajectories inside of region.
    :param disc: Discrepancy function to measure exactly on region.
    :return: red set of trajectories, blue set of trajectories, planted region, discrepancy function value.
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

    return red_in + red_out, blue_in + blue_out, plant_region


def plant_partial_halfplane(trajectories, r, p, q, eps):
    """
    Choose a random direction and expand the region along this direction till we contain r fraction of the total
    trajectory arc length.

    :param trajectories: List of list of points.
    :param r: double between 0 and 1
    :param p: double between 0 and 1
    :param q: double between 0 and 1
    :param eps: This defines the maximum difference between r and the fraction of the trajectories that are in the found region.
    :param disc: The discrepancy function to evaluate exactly on this region.
    :return: red set of trajectories, blue set of trajectories, the planted region, and the exact partial discrepancy.
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
    return red_in + red_out, blue_in + blue_out, plant_region


def plant_full_disk(trajectories, r, p, q):
    """
    Choose a point at random from a trajectory and then expand outward from there until we find a disk contains r fraction
    of the trajectories.

    :param trajectories: List of trajectories. This can be rearanged.
    :param r: Fraction of trajectories in region.
    :param p: Fraction of red trajectories outside of region.
    :param q: Fraction of red trajectories inside of region.
    :param disc: Discrepancy function to measure exactly on region.
    :return: red set, blue set, planted disk, maximum discrepancy.
    """
    traj = trajectories[random.randint(0, len(trajectories)-1)]
    origin = traj[random.randint(0, len(traj)-1)]

    trajectories.sort(key=lambda el: Trajectory(el).point_dist(origin))

    trajectories_ix = list(range(len(trajectories)))

    inside_disk = trajectories_ix[:int(r * len(trajectories))]
    red_in, blue_in = split_set(inside_disk, q)
    max_disk = Disk(origin[0], origin[1], Trajectory(trajectories[inside_disk[-1]]).point_dist(origin))
    del inside_disk

    outside_disk = trajectories_ix[int(r * len(trajectories)):]
    red_out, blue_out = split_set(outside_disk, p)
    del outside_disk

    red_ix = red_in + red_out
    blue_ix = blue_in + blue_out

    #sort the trajectories with respect to the ix.
    return [trajectories[i] for i in red_ix], [trajectories[i] for i in blue_ix], max_disk


def plant_partial_disk(trajectories, r, p, q, eps):
    """
    Choose a point at random from a trajectory and then expand outward from there. Computes the fraction of length inside the disk
    for each segment and then does bisection on this amount since it is a monotonic function to compute within eps fraction.

    :param trajectories: List of list of points.
    :param r: double between 0 and 1
    :param p: double between 0 and 1
    :param q: double between 0 and 1
    :param eps: This defines the maximum difference between r and the fraction of the trajectories that are in the found region.
    :param disc: The discrepancy function to evaluate exactly on this region.
    :return: red set of trajectories, blue set of trajectories, the planted region, and the exact partial discrepancy.
    """

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
    return red_in + red_out, blue_in + blue_out, max_disk


def distribution(points, p, scan_f, n, s, disc=DISC):
    while True:
        red, blue = split_set(points, p)

        net_set = my_sample(red, min(len(red), n)) + my_sample(blue, min(len(blue), n))
        m_sample = my_sample(red, min(len(red), s), red)
        b_sample = my_sample(blue, min(len(blue), s), blue)
        reg, val = scan_f(net_set, m_sample, b_sample, disc)

        yield val


def null_cdf(observations):
    values = sorted(observations)
    prob = [x / float(len(observations)) for x in range(len(observations))]
    return values, prob


def random_rect(points, r):
    """
    Plants a random rectangle containing r fraction of the points.

    :param points: List of points
    :param r: Fraction of points inside of planted region.
    :return: The planted region.
    """
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
                return Rectangle(ux, uy, lx, ly)
            el = ab.pop()
            inB.append(el)
            while len(ab) > 0 and sf(el) == sf(ab[-1]):
                el = ab.pop()
                inB.append(el)


def plant_rectangle(pts, r, p, q):
    """
    Create a set of red and blue points with a random rectangle planted containing r fraction of the points with
    q fraction of the points in the region being red and p fraction of the points outside of the region being red.

    :param pts: List of points.
    :param r: Fraction of points contained in the planted region
    :param p: Fraction of points outside region that red.
    :param q: Fraction of points inside region that are red.
    :return: red set, blue set, planted region.
    """

    rect = random_rect(pts, r)

    inside_rect = []
    outside_rect = []
    for pt in pts:
        if rect.contains(pt):
            inside_rect.append(pt)
        else:
            outside_rect.append(pt)

    red_in, blue_in = split_set(inside_rect, q)
    red_out, blue_out = split_set(outside_rect, p)
    return red_in + red_out, blue_in + blue_out, rect


def plant_disk(pts, r, p, q):
    """
    Create a set of red and blue points with a random disk planted containing r fraction of the points with
    q fraction of the points in the region being red and p fraction of the points outside of the region being red.

    :param pts: List of points.
    :param r: Fraction of points contained in the planted region
    :param p: Fraction of points outside region that red.
    :param q: Fraction of points inside region that are red.
    :return: red set, blue set, planted region.
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

    red_in, blue_in = split_set(inside_disk, q)
    red_out, blue_out = split_set(outside_disk, p)
    return red_in + red_out, blue_in + blue_out, max_disk


def plant_halfplane(pts, r, p, q):
    """
    Create a set of red and blue points with a random halfplane planted containing r fraction of the points with
    q fraction of the points in the region being red and p fraction of the points outside of the region being red.

    :param pts: List of points
    :param r: Fraction of points contained in the planted region
    :param p: Fraction of points outside region that red.
    :param q: Fraction of points inside region that are red.
    :return: red set, blue set, planted region.
    """

    def min_distance(pt, direc):
        return direc[0] * pt[0] + direc[1] * pt[1] 
    rand_direc = (random.gauss(0, 1), abs(random.gauss(0, 1)))

    ordered = sorted(pts, key=lambda x: min_distance(x, rand_direc))
    inside_halfplane = ordered[:int(r * len(ordered))]
    outside_halfplane = ordered[int(r * len(ordered)):]

    pt = ordered[int(r * len(ordered))]
    lc = pt[0] * rand_direc[0] + pt[1] * rand_direc[1]


    plant_region = Halfplane(Point(-rand_direc[0], -rand_direc[1], lc))

    red_in, blue_in = split_set(inside_halfplane, q)
    red_out, blue_out = split_set(outside_halfplane, p)
    return red_in + red_out, blue_in + blue_out, plant_region


def plant_partial_rectangle(trajectories, r, p, q, eps):
    """
    This plants a region containing r fraction of the total arc length of all trajectories. q fraction of the trajectories
    crossing this region are assigned to be anomalous and p fraction of trajectories not crossing this region are assigned
    to be anomalous.

    :param trajectories: List of list of points.
    :param r: double between 0 and 1
    :param p: double between 0 and 1
    :param q: double between 0 and 1
    :param eps: This defines the maximum difference between r and the fraction of the trajectories that are in the found region.
    :return: red set of trajectories, blue set of trajectories, the planted region.
    """
    trajectory_obj = [Trajectory(pts) for pts in trajectories]
    all_pts = uniform_sample(trajectory_obj, int(1 / eps ** 2 + 1), False)
    _, _, rect = plant_rectangle(all_pts, r, p, q)
    inside_rect = [traj for traj in trajectory_obj if rect.intersects_trajectory(traj)]
    outside_rect = [traj for traj in trajectory_obj if not rect.intersects_trajectory(traj)]
    red_in, blue_in = split_set(inside_rect, q)
    red_out, blue_out = split_set(outside_rect, p)

    return red_in + red_out, blue_in + blue_out, rect




def paired_plant_region(traj_start, traj_end, r, q, region_plant_f):
    """
    This plants a region where every trajectory
    completely outside or inside of the region has an endpoint chosen at random.
    Every trajectory with one endpoint inside the region has an endpoint chosen inside
    with probability q (exactly q fraction have one endpoint in the region)
    traj_start and traj_end should be the same length.

    :param traj_start: List of points.
    :param traj_end: List of points.
    :param r: Fraction of points in the region
    :param q: Fraction of points in the region that are anomalous
    :param region_plant_f: Scanning function to use to find the region (example max_disk)
    :return: Red planted set, blue planted set, and the planted region.
    """
    _, _, reg = region_plant_f(traj_start + traj_end, r, .5, q)

    flux_region = []
    out_region = []
    for st_pt, end_pt in zip(traj_start, traj_end):
        if reg.contains(st_pt) and not reg.contains(end_pt):
            flux_region.append((st_pt, end_pt))
        elif reg.contains(end_pt) and not reg.contains(st_pt):
            flux_region.append((end_pt, st_pt))
        else:
            if random.random() <= .5:
                out_region.append((st_pt, end_pt))
            else:
                out_region.append((end_pt, st_pt))


    q_fraction, remainder = split_set(flux_region, q)
    remainder = [(ep, sp) for (sp, ep) in remainder]

    red, blue = zip(*(q_fraction + remainder + out_region))
    return red, blue, reg



def plant_kernel_disk_region(pts, r, p, q, eps=.0001):
    """
    Plants a region using the bernoulli model with the given p and q parameters.
    Need the fraction in the q set to be r.

    Returns the bandwidth of the anomaly.
    :param pts:
    :param r:
    :param p:
    :param q:
    :param eps: The error on the planted r parameter.
    :return:
    """
    if not pts:
        return None
    ((min_x, min_y), (max_x, max_y)) = bounding_box(pts)


    def sq_dist(x, pt):
        return (x[0] - pt[0])**2 + (x[1] - pt[1])**2

    [p1, p2] = random.sample(pts, 2)
    alpha = random.random()
    px = p1[0] * alpha + p2[0] * (1 - alpha)
    py = p1[1] * alpha + p2[1] * (1 - alpha)

    #Find the correct bandwidth by doing bisection
    seeds = [random.random() for p in pts]
    bandwidth_upper = 2 * max((max_x - min_x), (max_y - min_y))
    bandwidth_lower = 0.0

    while True:
        bandwidth = (bandwidth_upper + bandwidth_lower) / 2
        def kernel(x):
            return math.exp(-sq_dist(x, (px, py)) / (2 * bandwidth**2))

        p_count = 0
        for i, pt in enumerate(pts):
            if seeds[i] < kernel(pt):
                p_count += 1

        if r - eps <= p_count / len(pts) <= r + eps:
            break
        if r + eps < p_count / len(pts):
            bandwidth_upper = bandwidth
        else:
            bandwidth_lower = bandwidth

    measured_set = []
    baseline_set = []
    p_count = 0
    q_count = 0
    for i, pt in enumerate(pts):
        if seeds[i] < kernel(pt):
            p_count += 1
            if random.random() < p:
                measured_set.append(pt)
            else:
                baseline_set.append(pt)
        else:
            q_count += 1
            if random.random() < q:
                measured_set.append(pt)
            else:
                baseline_set.append(pt)

    return measured_set, baseline_set, bandwidth, (px, py)



def disc_bernoulli_kern(measured, baseline, p, q, bandwidth, center):
    def sq_dist(x, pt):
        return (x[0] - pt[0])**2 + (x[1] - pt[1])**2
    def kernel(x):
        return math.exp(-sq_dist(x, center) / bandwidth**2)

    act_disc = 0
    for pt in measured:
        fr = kernel(pt)
        gr = p * fr + (1 - fr) * q
        act_disc += math.log(gr)
    for pt in baseline:
        fr = kernel(pt)
        gr = p * fr + (1 - fr) * q
        act_disc += math.log(1 - gr)

    scale = len(measured) / (len(measured) + len(baseline))
    null_val = len(measured) * math.log(scale) + len(baseline) * math.log(1 - scale)
    return act_disc - null_val



def close_region(region):
    """
    Ensures that this region starts and ends at the same point. This is useful for using trajectory methods
    for regions.
    :param region: List of points
    :return: list of points
    """
    if region:
        if Point(region[-1][0], region[-1][1], 1.0).approx_eq(region[0]):
            region.append(region[0])

def plant_full_disk_region(region_set, r, p, q):
    """
    Choose a point at random from a region and then expand a disk out from this point till this disk contains
    r fraction of all the regions.
    :param region_set: List of list of points. This argument is modified by adding a point at the end in some cases.
    :param r: double between 0 and 1
    :param p: double between 0 and 1
    :param q: double between 0 and 1
    :return: red set of regions, blue set of regions, the planted region
    """
    for reg in region_set:
        close_region(reg)
    return plant_full_disk(region_set, r, p, q)


try:
    from shapely.geometry import Polygon
    from shapely import geometry



    def measure_range_region(reg_geom, range_set, weights=None):
        weight_sum = 0
        if weights is None:
            weights = [1.0] * len(range_set)
        if isinstance(reg_geom, Disk):

            for i in range(len(range_set)):
                poly = Polygon(range_set[i])
                pt = reg_geom.get_origin()
                pt_tpl = geometry.Point(pt[0], pt[1])
                if poly.distance(pt_tpl) <= reg_geom.get_radius():
                    weight_sum += weights[i]
        elif isinstance(reg_geom, Halfplane):
            for i in range(len(range_set)):
                traj = Trajectory(range_set[i])
                if reg_geom.intersects_trajectory(traj):
                    weight_sum += weights[i]

        elif isinstance(reg_geom, Rectangle):
            rect = geometry.box(reg_geom.lowX(), reg_geom.lowY(), reg_geom.upX(), reg_geom.upY())
            for i in range(len(range_set)):
                poly = Polygon(range_set[i])
                if rect.intersects(poly):
                    weight_sum += weights[i]
        else:
            raise ValueError()
        return weight_sum
except:
    pass

try:
    import matplotlib.pyplot as plt
    import numpy as np

    def plot_points(ax, pts, c):
        xs = []
        ys = []
        for pt in pts:
            xs.append(pt[0] )
            ys.append(pt[1])
        ax.scatter(xs, ys, color=c, marker='.')

    def plot_kernel(ax, pts, pt, bandwidth, res=20, transform=None):
        (mnx, mny), (mxx, mxy) = bounding_box(pts)
        mxx = np.linspace(mnx, mxx, res)
        mxy = np.linspace(mny, mxy, res)
        xv, yv = np.meshgrid(mxx, mxy)

        def kernel(x, y):
            return np.exp(-(np.power(x - pt[0], 2.0) + np.power(y - pt[1], 2.0) ) / bandwidth**2)
        ax.contourf(xv, yv, kernel(xv, yv), alpha=.3, antialiased=True, transform=transform, cmap="Blues")
        ax.contour(xv, yv, kernel(xv, yv), linewidths=4.0, antialiased=True, transform=transform)

except:
    pass

def plant_full_square_region(regions, r, p, q, max_count=32):
    """
    Choose a point at random from a region and then expand a square out from this point till this region contains
    r fraction of all the regions.

    :param regions: List of list of points. These lists are modified. In some cases the regions will have a point appended at the end to close them.
    :param r: double between 0 and 1
    :param p: double between 0 and 1
    :param q: double between 0 and 1
    :param max_count: The maximum number of times we will attempt to find the right sized region.
    :return: red set of regions, blue set of regions, the planted region
    """
    for reg in regions:
        close_region(reg)
    return plant_full_square(regions, r, p, q, max_count)

def max_disk_trajectory(net, red_sample, blue_sample, min_disk_r, max_disk_r, alpha, disc, fast_disk=True):
    """
    Computes the highest discrepancy disk over a set of trajectories. Executes at multiple scales using the grid
    directional compression method and internally compresses the trajectories if fast_disk is enabled.

    :param net: A list of trajectories (scaled to be in a 0 by 1 box).
    :param red_sample: A list of trajectories (scaled to be in a 0 by 1 box).
    :param blue_sample: A list of trajectories (scaled to be in a 0 by 1 box).
    :param min_disk_r: The minimum disk radius to consider.
    :param max_disk_r: The maximum disk radius to consider.
    :param alpha: The spatial error with which to approximate the trajectories.
    :param disc: The discrepancy function to use.
    :param fast_disk: Default True.
    :return:
    """
    mx = -1
    curr_disk_r = max(min_disk_r, alpha)
    reg = None
    while True:

        chord_l = math.sqrt(4 * alpha * curr_disk_r - 2 * alpha * alpha)
        m_sample = [grid_direc_kernel(dp_compress(traj, alpha), chord_l, alpha) for traj in red_sample]
        b_sample = [grid_direc_kernel(dp_compress(traj, alpha), chord_l, alpha) for traj in blue_sample]
        pt_net = [grid_direc_kernel(dp_compress(traj, alpha), chord_l, alpha) for traj in net]
        m_sample = list(trajectories_to_labels(m_sample))
        b_sample = list(trajectories_to_labels(b_sample))
        net_set = list(trajectories_to_labels(pt_net))


        new_reg, new_mx = max_disk_scale_labeled(net_set, m_sample, b_sample, fast_disk, curr_disk_r, disc)
        if new_mx > mx:
            reg = new_reg
            mx = new_mx
        curr_disk_r *= 2
        if curr_disk_r >= max_disk_r:
            break
    return reg, mx


def max_disk_trajectory_fixed(net, m_sample, b_sample, min_disk_r, max_disk_r,  disc, fast_disk=True):
    """
    Computes the highest discrepancy disk over a set of trajectories. Executes at multiple scales, but uses whatever set
    of points the trajectories have been compressed with.

    :param net: A list of trajectories (scaled to be in a 0 by 1 box).
    :param m_sample: A list of trajectories (scaled to be in a 0 by 1 box).
    :param b_sample: A list of trajectories (scaled to be in a 0 by 1 box).
    :param min_disk_r: The minimum disk radius to consider.
    :param max_disk_r: The maximum disk radius to consider.
    :param disc: The discrepancy function to use.
    :param fast_disk: Default True.
    :return: A tuple of the maximum disk and the corresponding maximum value.
    """
    mx = -1
    curr_disk_r = min_disk_r
    reg = None
    while curr_disk_r < max_disk_r:
        new_reg, new_mx = max_disk_scale_labeled(net, m_sample, b_sample, fast_disk, curr_disk_r, disc)
        if new_mx > mx:
            reg = new_reg
            mx = new_mx
        curr_disk_r *= 2
    return reg, mx


def max_disk_region(net, red_sample, red_weight, blue_sample, blue_weight, min_disk_r, max_disk_r, alpha, disc, fast_disk=True):
    """
    Computes the highest discrepancy disk over a set of trajectories. Executes at multiple scales using the grid
    directional compression method and internally compresses the trajectories if fast_disk is enabled.

    :param net: A list of trajectories
    :param red_sample: A list of trajectories
    :param blue_sample: A list of trajectories
    :param min_disk_r: The minimum disk radius to consider.
    :param max_disk_r: The maximum disk radius to consider.
    :param alpha: The spatial error with which to approximate the trajectories.
    :param disc: The discrepancy function to use.
    :param fast_disk: Default True.
    :return:
    """
    mx = -1
    curr_disk_r = max(min_disk_r, alpha)
    reg = None
    while True:

        m_sample = [polygon_grid_hull(reg, alpha, curr_disk_r) for reg in red_sample]
        b_sample = [polygon_grid_hull(reg, alpha, curr_disk_r) for reg in  blue_sample]
        pt_net = [polygon_grid_hull(reg, alpha, curr_disk_r) for reg in net]

        m_sample = [[LPoint(i, w, pt[0], pt[1], 1.0) for pt in reg] for reg, w, i in zip(m_sample, red_weight, range(len(m_sample)))]
        b_sample = [[LPoint(i, w, pt[0], pt[1], 1.0) for pt in reg] for reg, w, i in zip(b_sample, blue_weight, range(len(m_sample)))]
        m_sample = list(itertools.chain.from_iterable(m_sample))
        b_sample = list(itertools.chain.from_iterable(b_sample))
        net_set = list(trajectories_to_labels(pt_net))
        new_reg, new_mx = max_disk_scale_labeled(net_set, m_sample, b_sample, fast_disk, curr_disk_r, disc)
        if new_mx > mx:
            reg = new_reg
            mx = new_mx
        curr_disk_r *= 2
        if curr_disk_r >= max_disk_r:
            break
    return reg, mx


