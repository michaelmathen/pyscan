
.. _utility-label:

Utility
===================

These methods are mostly utility methods written in python in pyscan.py.

.. automodule:: pyscan
     :members:


.. py:function:: evaluate_halfplane(reg, mpts, bpts, disc_f)

.. py:function:: evaluate_halfplane_labeled(reg, lmpts, lbpts, disc_f)

.. py:function:: evaluate_halfplane_trajectory(reg, mtrajs, btrajs, disc_f)

.. py:function:: evaluate_disk(reg, mpts, bpts, disc_f)

.. py:function:: evaluate_disk_alt(reg, mpts, bpts, disc_f)

.. py:function:: evaluate_disk_labeled(reg, lmpts, lbpts, disc_f)

.. py:function:: evaluate_disk_trajectory(reg, mtrajs, btrajs, disc_f)

.. py:function:: evaluate_rectangle(reg, mpts, bpts, disc_f)

.. py:function:: evaluate_rectangle_labeled(reg, lmpts, lbpts, disc_f)

.. py:function:: evaluate_rectangle_trajectory(reg, mtrajs, btrajs, disc_f)

   These functions can be used to evaluate different regions on different kinds of data sets.

.. py:function:: size_region(fraction)
  
   Creates a discrepancy function that can be used to find a region containing a certain fraction of the points. 
   
   :param fraction: A double between 0 and 1.
   :rtype: A discrepancy function.

.. py:function:: evaluate(disc_f, m, m_total, b, b_total)
  
   A utility function for explicitly evaluating a discrepancy function object.
   
   :param disc_f: A double between 0 and 1.
   :param m: A double between 0 and m_total
   :param m_total: A double.
   :param b: A double between 0 and b_total.
   :param b_total: A double.
   :rtype: A double returned by the discrepancy function.
