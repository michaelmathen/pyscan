
Api
====================

The base data types.

* :any:`Point`
* :any:`WPoint`
* :any:`LPoint`
* :any:`Point3`
* :any:`WPoint3`
* :any:`LPoint3`
* :any:`Trajectory`

There are several types of regions of interest.

* :any:`Halfplane`
* :any:`Halfspace`
* :any:`Disk`
* :any:`Rectangle`

We have a variety of types of scanning algorithms. Unlabeled scanning algorithms solve the classic bichromatic discrepancy problem or the approximate scan statistic problem. The labeled versions allow the user to mark points depending on the object they came from. This allows region scanning or trajectory scanning problems to be reinterpreted as labeled scanning problems. We also have 
recently begun work on kernelized scan statistics. In this project formulation the region is defined as a kernel with a certain bandwidth and the scan statistic must be redefined to deal with correlation that can occur between points. This last method is under heavy active development, so expect changes.

* :ref:`scanning-label`
* :ref:`labeled-scanning-label`
* :ref:`smooth-scanning-label`

We also have some older codes that were used for experiments in :cite:`AMPVZ06` that have been wrapped for use in the python framework.

* :ref:`jeff-codes`

Utility Methods :ref:`utility-label`.












