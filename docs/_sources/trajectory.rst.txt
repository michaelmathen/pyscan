Trajectory
===================

.. py:currentmodule:: pyscan

.. py:class:: pyscan.Trajectory(pts)
   
   This object represents a trajectory.

   :param pts: List containing Point objects

.. py:method:: Trajectory.point_dist(pt)
   
   Computes the minimum distance from this point to the trajectory.

   :param pt: Point
   :rtype: double

.. py:method:: Trajectory.get_weight()
   
   Computes the weight of the trajectory. By default this will be 1 unless you create a weighted trajectory.

   :rtype: double

.. py:method:: Trajectory.get_length()
   
   Computes the weight of the trajectory. By default this will be 1 unless you create a weighted trajectory.

   :rtype: double

.. py:method:: Trajectory.get_pts()
   
   Returns the trajectories internal point list.

   :rtype: List containing Point objects

.. py:class:: pyscan.WTrajectory(weight, pts)
   
   This object represents a weighted trajectory.

   :param weight: double
   :param pts: List containing Point objects


Trajectory Simplification
============================
These methods can be used to take get spatial approximations for trajectories to allow for the full scanning model to be applied.

.. py:function:: dp_compress(traj, alpha)
   
   Douglas-Peuker trajectory simplification method :cite:`RAMER72`.

   :param traj: List of points.
   :param alpha: The pruning distance in the DP algorithm.
   :rtype: A list of points corresponding to the simplified trajectory.


.. py:function:: grid_kernel(traj, alpha)
   
   Grids the trajectory using a :math:`1/\alpha \times 1/\alpha` grid and picks the center point of each cell the trajectory crosses. 

   :param traj: List of points.
   :param alpha: double
   :rtype: A list of points corresponding to the simplified trajectory.

.. py:function:: grid_trajectory(traj, alpha)
   
   Grids the trajectory using a :math:`1/\alpha \times 1/\alpha` grid and picks the points where the trajectory crosses the boundaries of the cell. 

   :param traj: List of points.
   :param alpha: double
   :rtype: A list of points corresponding to the simplified trajectory.

.. py:function:: grid_direc_kernel(traj, r, alpha)
   
   Grids the trajectory using a :math:`1/r \times 1/r` grid and then applies a directional kernel to each cell with :math:`\alpha` error. This simplification get error
   of :math:`\alpha` for disks with radii greater than :math:`r`.

   :param traj: List of points.
   :param r: double
   :param alpha: double
   :rtype: A list of points corresponding to the simplified trajectory.

.. py:function:: halfplane_kernel(traj, alpha)
   
   Applies a directional kernel with :math:`\alpha` error. This simplification get error of :math:`\alpha` for halfplanes.

   :param traj: List of points.
   :param alpha: double
   :rtype: A list of points corresponding to the simplified trajectory.

.. py:function:: hull(traj)
   
   Takes the convex hull of the trajectory. This simplification has zero error for halfplanes and can significantly speed up scanning for labeled data.

   :param traj: List of points.
   :rtype: A list of points corresponding to the simplified trajectory.

.. py:function:: lifting_kernel(traj, alpha)
   
   Takes the convex hull of the trajectory lifted to a 3d paraboloid. This simplification has :math:`alpha` error for disks, but currently has stability issues.

   :param traj: List of points.
   :param alpha: double
   :rtype: A list of points corresponding to the simplified trajectory.


.. py:function:: uniform_sample_error(traj, alpha, endpoints)
   
   Chooses uniformly randomly :math:`L/\alpha` points on the trajectory where :math:`L` is the arc length. If the endpoint variable is true then this method 
   also takes the endpoints of the trajectory.

   :param traj: List of points.
   :param alpha: double
   :param endpoints: boolean
   :rtype: A list of points corresponding to the simplified trajectory. 

.. py:function:: even_sample_error(traj, alpha, endpoints)
   
   Chooses a point ever :math:`\alpha` distance apart. If the endpoint variable is true then this method 
   also takes the endpoints of the trajectory.

   :param traj: List of points.
   :param alpha: double
   :param endpoints: boolean
   :rtype: A list of points corresponding to the simplified trajectory.

Trajectories to Points
============================
These algorithms can be used to convert the trajectory to a point set that can then be fed into one of the approximate scanning methods.

.. py:function:: block_sample(trajectories, s, endpoints)
   
   Chooses s points such that one point is chosen at random inside of each arc length segment of length :math:`L/s` where L is the summed length of all the 
   trajectories in trajectories. If the endpoint variable is true then this method 
   also takes the endpoints of the trajectory.

   :param trajectories: A list of Trajectories
   :param alpha: integer
   :param endpoints: boolean
   :rtype: A list of weighted points

.. py:function:: uniform_sample(trajectories, s, endpoints)
   
   Uniformly samples :math:`s` points on the trajectories. If the endpoint variable is true then this method 
   also takes the endpoints of the trajectory.

   :param traj: A list of Trajectories
   :param alpha: integer
   :param endpoints: boolean
   :rtype: A list of weighted points

.. py:function:: even_sample(trajectories, s, endpoints)
   
   Chooses :math:`s` points such that one point is every arc length segment of length :math:`L/s` where L is the summed length of all the 
   trajectories in trajectories. If the endpoint variable is true then this method 
   also takes the endpoints of the trajectory.

   :param traj: A list of Trajectories
   :param alpha: integer
   :param endpoints: boolean
   :rtype: A list of weighted points
 

