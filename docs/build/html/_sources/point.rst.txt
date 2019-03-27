
Point
====================
.. py:currentmodule:: pyscan

.. py:class:: pyscan.Point(a, b, c)
   
   This is a point class. Internally the data is represented using a homogenous coordinate system. For those not familiar this means that a homogenous coordinate :math:`(a, b, c)` in standard coordinates would be :math:`(a / c, b / c)`). These points can also be used to represent lines in which case the homogenous coordinate :math:`(a, b, c)` is the line :math:`a x + b y + c = 0`.
  
   :param a: double
   :param b: double
   :param c: double

>>> import pyscan
>>> p1 = pyscan.Point(1.0, 2.0, 1.0)
>>> p2 = pyscan.Point(5.0, 10.0, 1.0) 
>>> p1.approx_eq(p2)
False
>>> p1.parallel_lte(p2)
True
>>> p1.evaluate(p2)
26.0
>>> p1[0]
1.0


.. py:method:: Point.approx_eq(p_other)
   
   Checks to see if these points are approximately the same. Internally compares floats using float distance.

   :param p_other: Another Point object.
   :rtype: True or False.

.. py:method:: Point.__getitem__(index)

   Returns the value at the standard coordinate index.

   :param index: An integer.
   :rtype: A double.

.. py:method:: Point.get_coord(index)
   
   Returns the homogenous coordinate at index.

   :param index: An integer.
   :rtype: A double.

.. py:method:: Point.above(pt)

   Returns whether this point is above the line represented by pt. Note that a line and point can have an upward or downward orientation.

   :param pt: A line (or point).
   :rtype: boolean.

.. py:method:: Point.above_closed(pt)

   Returns whether this point is above the line or contained inside of the line represented by pt. Note that a line and point can have an upward or downward orientation.

   :param pt: A line (or point).
   :rtype: boolean.

.. py:method:: Point.below(pt)

   Returns whether this point is below the line represented by pt. Note that a line and point can have an upward or downward orientation.

   :param pt: A line (or point).
   :rtype: boolean.

.. py:method:: Point.below_closed(pt)

   Returns whether this point is below the line or contained inside of the line represented by pt. Note that a line and point can have an upward or downward orientation.

   :param pt: A line (or point).
   :rtype: boolean.

.. py:method:: Point.crosses(p1, p2)

   Returns whether the line represented by the homogenous coordinates in point crosses the line segment between p1 and p2. In the dual formulation this can be thought of as 
   checking if this point lies inside of the double wedge defined by the lines p1 and p2.

   :param p1: A point(or line).
   :param p2: A point(or line).
   :rtype: boolean.

.. py:method:: Point.evaluate(pt)

   Takes the dot product between this point's and the argument's homogenous coordinates.
 
   :param p1: A point(or line).
   :rtype: double

.. py:method:: Point.orient_down(ix)

   Ensures that the line's normal is oriented such that ix coordinate of the normal is negative.
 
   :param ix: Integer between 0 and 2.
   :rtype: Point

.. py:method:: Point.orient_up(ix)

   Ensures that the line's normal is oriented such that ix coordinate of the normal is positive.
 
   :param ix: Integer between 0 and 2.
   :rtype: Point

.. py:method:: Point.dist(pt)

   Returns the euclidean distance between this point and pt.
 
   :param pt: Point object
   :rtype: double

.. py:method:: Point.pdot(pt)

   Dot product using the normalized coordinates. 
 
   :param pt: Point object
   :rtype: double

.. py:method:: Point.parallel_lt(pt)

   Checks to see if this line is parallel to the line represented by pt and also happens to be less than it in terms of their orientations.
 
   :param pt: Point object
   :rtype: bool

.. py:method:: Point.parallel_lte(pt)

   Checks to see if this line is parallel to the line represented by pt and also happens to be less than or equal to it in terms of their orientations.
 
   :param pt: Point object
   :rtype: bool


.. py:class:: pyscan.WPoint(weight, a, b, c)

   This inherits all the methods that :any:`Point` has, but also has positive weight.

   :param weight: double
   :param a: double
   :param b: double
   :param c: double

.. py:method:: WPoint.get_weight()

   Returns this points weight.
 
   :rtype: double

.. py:class:: pyscan.LPoint(label, weight, a, b, c)

   This inherits all the methods that :any:`WPoint` has, but also has a integer label.

   :param label: A non negative integer.
   :param weight: double
   :param a: double
   :param b: double
   :param c: double

.. py:method:: LPoint.get_label()

   Returns this points label.
 
   :rtype: integer


.. py:method:: LPoint.get_label()

   Returns this points label.
 
   :rtype: integer

.. py:class:: pyscan.Point3(a, b, c, d)

   Same methods as :any:`Point`, but for 3d points.

.. py:class:: pyscan.WPoint3(weight, a, b, c, d)

.. py:class:: pyscan.LPoint3(label, weight, a, b, c, d)

   
