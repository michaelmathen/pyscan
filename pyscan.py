from libpyscan import *

def trajectories_to_flux(trajectories):
    return [t[0] for t in trajectories], [t[1] for t in trajectories]

