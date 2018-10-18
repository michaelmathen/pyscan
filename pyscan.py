from libpyscan import *
import random


def trajectories_to_flux(trajectories):
    return [t[0] for t in trajectories], [t[1] for t in trajectories]


def evaluate_range(range, mp, bp, scan_f):
    if type(range) is Disk:
        return evaluate_disk(range, mp, bp, scan_f)
    elif type(range) is Halfplane:
        return evaluate_halfspace(range, mp, bp, scan_f)
    else:
        raise "Not a valid option"



def split_set(pts, rate):
    red_set = random.sample(pts, len(pts) * rate)
    red_set_set = set(red_set)
    return red_set, [item for item in pts if item not in red_set_set]


"""
This plants a region where every trajectory:
Completely outside or inside of the region has an endpoint chosen at random.
Every trajectory with one endpoint inside the region has an endpoint chosen inside
with probability q (exactly q fraction have one endpoint in the region)

r controls how many points the region contains.

"""
def paired_plant_region(traj_start, traj_end, r, q, eps, scan_f):

    while True:

        net_size = int( 1 / (eps * r) + 1)
        s_size = int(1 / (eps * r) ** 2 + 1)

        net = random.sample(traj_start + traj_end, net_size)
        sample = random.sample(traj_start + traj_end, s_size)
        sample = [WPoint(1.0, p[0], p[1], 1.0) for p in sample]
        disc = pyscan.size_region(r)
        reg, _ = scan_f(net, sample, [], disc)

        val = evaluate_range(reg, traj_start + traj_end, [], disc)
        if 1 - val < eps:
            break

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

    return zip(*(q_fraction + q_fraction_o + remainder + remainder_o))


"""
This takes a scanning function and two point sets and then computes a planted region that contains some fraction r
of the points with some tolerance.
This then computes a red and blue set of points based on the planted region.

Inside the region q fraction of points are red.
Outside the region p fraction of points are red
"""
def plant_region(points, r, p, q, eps, scan_f):

    while True:

        net_size = int(1 / (eps * r) + 1)
        s_size = int(1 / (eps * r) ** 2 + 1)

        net = random.sample(points, net_size)
        sample = random.sample(points, s_size)
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
    return red_in + red_out, blue_in + blue_out


def distribution(points, p, scan_f, n, s, disc=pyscan.DISC):
    for i in range(count):
        red, blue = split_set(pts, p)

        net_set = random.sample(min(len(red), n), red) + random.sample(min(len(blue), n), blue)
        m_sample = random.sample(min(len(red), s), red)
        b_sample = random.sample(min(len(blue), s), blue)
        reg, val = scan_f(net_set, m_sample, b_sample, disc)

        yield val

