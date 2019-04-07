Ranges
===============================
This is the api for all the range objects.

.. py:class:: pyscan.Halfplane(pt)

   This is a range representing a halfplane. Internally it uses a homogenous coordinate and is usually oriented downwards.

   :param pt: Point

.. py:method:: Halfplane.get_coords()
   
   Returns the internal :any:`Point` representation of this :any:`Halfplane`.

   :rtype: Point

.. py:method:: Halfplane.contains(pt)

   Checks to see if the halfplane contains the point. Internally uses :any:`Point.below_closed`.

   :param pt: Point
   :rtype: boolean

.. py:method:: Halfplane.intersects_segment(p1, p2)

   Checks to see if the halfplane intersects the line segment represented by p1 and p2. Internally uses :any:`Point.crosses`.

   :param p1: Point
   :param p2: Point
   :rtype: boolean

.. py:method:: Halfplane.intersects_trajectory(traj)

   Checks to see if the halfplane intersects the trajectory represented by traj.

   :param traj: Trajectory
   :rtype: boolean

.. py:class:: pyscan.Halfspace(pt)

   This is a range representing a halfspace in 3d. Internally it uses a homogenous coordinate and is usually oriented downwards.

   :param pt: Point3

.. py:method:: Halfspace.get_coords()
   
   Returns the internal :any:`Point3` representation of this :any:`Halfplane`.

   :rtype: Point3

.. py:method:: Halfplane.contains(pt)

   Checks to see if the halfplane contains the point. Internally uses :any:`Point.below_closed`.

   :param pt: Point3
   :rtype: boolean

.. py:method:: Halfspace.intersects_segment(p1, p2)

   Checks to see if the halfplane intersects the line segment represented by p1 and p2. Internally uses :any:`Point.crosses`.

   :param p1: Point3
   :param p2: Point3
   :rtype: boolean


.. py:class:: pyscan.Disk(origin, radius)

   This is a range representing a disk. Internally it uses a point for its center and a radius around that center.

   :param origin: Point
   :param radius: double

.. py:method:: Disk.getOrigin()
   
   Returns the internal :any:`Point` representation of the disk center.

   :rtype: Point

.. py:method:: Disk.getRadius()
   
   Returns the radius of this disk.

   :rtype: double

.. py:method:: Disk.contains(pt)

   Checks to see if the Disk contains the point.

   :param pt: Point
   :rtype: boolean

.. py:method:: Disk.intersects_segment(p1, p2)

   Checks to see if the Disk intersects the line segment represented by p1 and p2.

   :param p1: Point
   :param p2: Point
   :rtype: boolean

.. py:method:: Disk.intersects_trajectory(traj)

   Checks to see if the disk intersects the trajectory represented by traj.

   :param traj: Trajectory
   :rtype: boolean


.. py:class:: pyscan.Rectangle(upX, upY, lowX, lowY)

   This is a range representing a rectangle. (upX, upY) is the upper right corner and (lowX, lowY) is the lower left corner.

   :param upX: double
   :param upY: double
   :param lowX: double
   :param lowY: double


.. py:method:: Rectangle.lowX()

   :rtype: double

.. py:method:: Rectangle.lowY()

   :rtype: double

.. py:method:: Rectangle.upX()

   :rtype: double

.. py:method:: Rectangle.upY()

   :rtype: double

.. py:method:: Rectangle.contains(pt)

   Checks to see if the Rectangle contains the point.

   :param pt: Point
   :rtype: boolean

.. py:method:: Rectangle.intersects_segment(p1, p2)

   Checks to see if the Rectangle intersects the line segment represented by p1 and p2.

   :param p1: Point
   :param p2: Point
   :rtype: boolean

.. py:method:: Rectangle.intersects_trajectory(traj)

   Checks to see if the rectangle intersects the trajectory represented by traj.

   :param traj: Trajectory
   :rtype: boolean

