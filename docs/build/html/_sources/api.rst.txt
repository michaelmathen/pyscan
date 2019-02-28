

Python Functions 
===================

.. automodule:: pyscan
     :members:

.. py:currentmodule:: pyscan

C++ Functions
====================

.. py:function:: pyscan.max_halfplane(net, mpts, bpts, disc_f)

   Scans the set of halfplanes defined by the points in the net and maximizes the function disc_f with respect to the fraction of points contained in mpts and bpts. Function takes :math:`O(ns\log n)` where n is the net size and s is the size of the mpts and bts.

   :param net: List of Points
   :param mpts: List of WPoints
   :param bpts: List of WPoints
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found halfplane and the maximum disc_f value

.. py:function:: pyscan.max_halfplane_fast(r, mpts, bpts, disc_f)

   HERE BE DRAGONS!!! This method is still experimental and somewhat broken. This method uses a different method to scan the set of potential halfplanes with respect to the fraction of points contained in mpts and bpts. This code could potentially be much faster than the standard max_halfplane method.

   :param r: Number of regions asymptopically to consider.
   :param mpts: List of WPoints
   :param bpts: List of WPoints
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found halfplane and the maximum disc_f value

.. py:function:: pyscan.max_halfspace(net, mpts, bpts, disc_f)

   Scans the set of halfspaces defined by the points in the net and maximizes the function disc_f with respect to the fraction of points contained in mpts and bpts. Function takes :math:`O(n^2s\log n)` where n is the net size and s is the size of the mpts and bts.


   :param net: List of Point3s
   :param mpts: List of WPoint3s
   :param bpts: List of WPoint3s
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found halfspace and the maximum disc_f value

.. py:function:: pyscan.max_halfplane_labeled(net, lmpts, lbpts, disc_f)
  
   Scans the set of halfplanes defined by the points in the net and maximizes the function disc_f with respect to the labeled sets of points lmpts and lbpts. If two points with the same label are in the same region then the two points only contribute one of their weights to the region. This algorithm runs in time :math:`O(ns \log n)` where n is the net size and s is the size of the lmpts and lbpts.

   :param net: List of Points
   :param mpts: List of LPoints
   :param bpts: List of LPoints
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found halfplane and the maximum disc_f value

.. py:function:: pyscan.max_halfspace_labeled(net, lmpts, lbpts, disc_f)
  
   Scans the set of halfspaces defined by the points in the net and maximizes the function disc_f with respect to the labeled sets of points lmpts and lbpts. If two points with the same label are in the same region then the two points only contribute one of their weights to the region. This algorithm runs in time :math:`O(n^2s \log n)` where n is the net size and s is the size of the lmpts and lbpts.

   :param net: List of Point3s
   :param mpts: List of LPoint3s
   :param bpts: List of LPoint3s
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found halfspace and the maximum disc_f value

.. py:function:: pyscan.ham_tree_sample(pts, size)

   Takes the list of pts and computes a ham tree sample :cite:`MP18a`. This method can shrink the size of the sample significantly while preserving the error with respect to halfplanes. This serves as a useful preprocessing step for speeding up halfplane scanning. In theory sample sizes go from :math:`O(\frac{1}{\varepsilon^2}` to :math:`O(\frac{1}{\varepsilon^{1.53..}}\log^{.766}\frac{1}{\varepsilon})` where :math:`\varepsilon` is the absolute fraction of misplaced points in a halfplane.

   :param pts: List of WPoints
   :rtype: List of WPoints

.. py:function:: pyscan.max_disk(net, mpts, bpts, disc_f)

   Scans the set of disks defined by the points in the net and maximizes the function disc_f with respect to the fraction of points contained in mpts and bpts. Function takes :math:`O(n^2s\log n)` where n is the net size and s is the size of the mpts and bts.


   :param net: List of Points
   :param mpts: List of WPoints
   :param bpts: List of WPoints
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found disk and the maximum disc_f value

.. py:function:: pyscan.max_disk_scale(net, mpts, bpts, min_res, disc_f)

   Scans the set of disks defined by the points in the net, between the min_res and 2 * min_res and maximizes the function disc_f with respect to the fraction of points contained in mpts and bpts. Internally we use the disk radius restriction to ignore far away points and prune out many potential disks. This can significantly decrease the amount of time the method takes to find the best region. In worst case this method will still take :math:`O(n^2s\log n)` where n is the net size and s is the size of the mpts and bts, but usually runtime will be much less.


   :param net: List of Points
   :param mpts: List of WPoints
   :param bpts: List of WPoints
   :param min_res: A float that corresponds to minimum radius to consider.
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found disk and the maximum disc_f value


.. py:function:: pyscan.max_disk_labeled(net, lmpts, lbpts, disc_f)

   Scans the set of disks defined by the points in the net and maximizes the function disc_f with respect to the fraction of points contained in lmpts and lbpts. Points with identical labels are not double counted if they are contained in the same region. In worst case this method will still take :math:`O(n^2s\log n)` where n is the net size and s is the size of the lmpts and lbts.


   :param net: List of Points
   :param lmpts: List of LPoints
   :param lbpts: List of LPoints
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found disk and the maximum disc_f value

.. py:function:: pyscan.max_disk_scale_labeled(net, lmpts, lbpts, compress, min_res, disc_f)

   Scans the set of disks defined by the points in the net and maximizes the function disc_f with respect to the fraction of points contained in lmpts and lbpts. Only considers disks with radii between min_res and 2 * min_res. This method internally maps the disk scanning problem to halfplane scanning. Using this mapping we can compress sets of points with the same labels. This compression step can significantly improve the runtime and doesn't incurr any extra error. In worst case this method will still take :math:`O(n^2s\log n)` where n is the net size and s is the size of the lmpts and lbts, but in practice this method should be much faster than this.


   :param net: List of Points
   :param lmpts: List of LPoints
   :param lbpts: List of LPoints
   :param compress: Turns the compression step on or off.
   :param min_res: Defines the radius range to consider.
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found disk and the maximum disc_f value

.. py:function:: pyscan.max_disk_scale_labeled_alt(net, lmpts, lbpts, min_res, disc_f)

   Scans the set of disks defined by the points in the net and maximizes the function disc_f with respect to the fraction of points contained in lmpts and lbpts. Only considers disks with radii between min_res and 2 * min_res. This method internally directly scans disks and lacks the compression step found in the previous algorihtm, so it is slighlty slower in practice then the previous method with compression turned off.  In worst case this method will still take :math:`O(n^2s\log n)` where n is the net size and s is the size of the lmpts and lbts, but in practice this method should be much faster than this.


   :param net: List of Points
   :param lmpts: List of LPoints
   :param lbpts: List of LPoints
   :param min_res: Defines the radius range to consider.
   :param disc_f: Discrepancy function to maximize over.
   :rtype: Tuple of the found disk and the maximum disc_f value

.. py:function:: pyscan.max_subgrid(grid, disc_f)
   
   Finds the maximum subgrid exactly by enumerating all the subgrids and computing disc_f on each subgrid. Takes :math:`O(n^4)` where the grid is :math:`n \times n`.
    
   :param grid: A Grid object
   :param disc_f: Discrepancy function to maximize over.
   :rtype: The maximum subgrid found.

.. py:function:: pyscan.max_subgrid_linear(grid, x, y)

   Finds the maximum subgrid with respect to the linear function :math:`f(m, b) = xm + yb` exactly by using the Kadane algorithm. Takes :math:`O(n^3)` where the grid is :math:`n \times n`.

   :param grid: A Grid object
   :param x: Double defines linear function.
   :param y: Double defines linear function.
   :rtype: The maximum subgrid found.

.. py:function:: pyscan.max_subgrid_linear_theory(grid, r, x, y)

   Finds an approximate maximum subgrid with respect to the linear function :math:`f(m, b) = xm + yb` using the algorithm from :cite:`MP18b`. Takes :math:`O(n^2 + r^2 \log r)` where the grid is :math:`n \times n`.
   :param grid: A Grid object
   :param r: Determines how frequently to approximate the subgrids. See :cite:`MP18b` for details.
   :param x: Double defines linear function.
   :param y: Double defines linear function.
   :rtype: The maximum subgrid found.

.. py:function:: pyscan.max_subgrid_convex(grid, eps, disc_f) 

   Approximately finds the maximum subgrid of the grid with error eps by successively computing the maximum subgrid with respect to various linear function using the Kadane algorithm. Takes :math:`O(\frac{1}{\sqrt{\varepsilon}} n^3)` where the grid is :math:`n \times n`.

   :param grid: A Grid object
   :param eps: The absolute error to incur.
   :param disc_f: Discrepancy function to maximize over.
   :rtype: The maximum subgrid found.

.. py:function:: pyscan.max_subgrid_convex_theory(grid, eps, disc_f)

   Approximately finds the maximum subgrid of the grid with error eps by successively computing the maximum subgrid with respect to various linear function using the algorithm from :cite:`MP18b`. Takes :math:`O(\frac{1}{\sqrt{\varepsilon}} n^2 \log n)` where the grid is :math:`n \times n`. This method is depreciated and will be replaced shortly.

   :param grid: A Grid object
   :param eps: The absolute error to incur.
   :param disc_f: Discrepancy function to maximize over.
   :rtype: The maximum subgrid found.

.. py:function:: pyscan.max_rectangle(mpts, bpts, eps, x, y)

.. py:function:: pyscan.max_rectangle_heap(mpts, bpts, x, y)

   Implements the algorithm from :cite:`APV06,DGM95` to scan the data set in :math:`s^2 \log s` time. We use a Treap based rotation scheme to keep the tree balanced.

   :param mpts: List of WPoints.
   :param bpts: List of WPoints.
   :param x: Linear parameter.
   :param y: Linear parameter.
   :rtype: The maximum rectangle found.

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

.. py:function:: dp_compress(traj, eps)
   
   Douglas-Peuker trajectory simplification method :cite:`RAMER72`.

   :param traj: The trajectory to simplify.
   :param eps: The pruning distance in the DP algorithm.
   :rtype: A list of points corresponding to the simplified trajectory.
 py::def("dp_compress", &pyscan::dp_compress);
    //This grids the trajectory and assigns a single point to each cell.
    py::def("grid_kernel", &pyscan::approx_traj_grid);
    py::def("grid_trajectory", &pyscan::grid_traj);
    //This grids the trajectory and creates an alpha hull in each one.
    py::def("grid_direc_kernel", &pyscan::approx_traj_kernel_grid);
    //This is for 2d eps-kernel useful for halfspaces.
    py::def("halfplane_kernel", pyscan::approx_hull);
    py::def("convex_hull", pyscan::graham_march);
    //This is a 3d eps-kernel for disks.
    py::def("lifting_kernel", &pyscan::lifting_coreset);

    py::def("coreset_error_halfplane", &pyscan::error_halfplane_coreset);
    py::def("coreset_error_disk", &pyscan::error_disk_coreset);
    py::def("hull", pyscan::graham_march);




