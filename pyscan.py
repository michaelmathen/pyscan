from libpyscan import *

def trajectories_to_flux(trajectories):
    return [t[0] for t in trajectories], [t[1] for t in trajectories]


def evaluate_range(range, mp, bp, scan_f):
    if type(range) is libpyscan.Disk:
        return libpyscan.evaluate_disk(range, mp, bp, scan_f)
    elif type(range) is libpyscan.Halfplane:
        return libpyscan.evaluate_halfspace(range, mp, bp, scan_f)
    else:
        raise "Not a valid option"