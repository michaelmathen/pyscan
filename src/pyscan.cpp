//
// Created by mmath on 7/7/17.
//

#include <functional>
#include <iostream>
#include <typeinfo>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "KernelScanning.hpp"
#include "JeffCodes.hpp"

#include "ConvexHull.hpp"
#include "Segment.hpp"
#include "RectangleScan.hpp"
#include "HalfSpaceScan.hpp"
#include "DiskScan.hpp"
#include "FunctionApprox.hpp"
#include "TrajectoryScan.hpp"
#include "ConvexHull.hpp"
#include "TrajectoryCoreSet.hpp"
#include "PartitionSample.hpp"
#include "RegionCoreSet.hpp"
#include "SatScan.hpp"


#define PY_WRAP(FNAME) py::def("FNAME", &pyscan:: FNAME)

namespace py = pybind11;


namespace pyscan {

    pyscan::Segment to_Segment(Point<2> const& pt) {
        return Segment(HalfSpace<2>(pt));
    }

    auto sized_region(double size) -> std::function<double(double, double, double, double)> {
        /*
         * Useful for finding a region of a certain size.
         */
        return [size] (double m, double m_total, double b, double b_total) {
            (void)b;
            (void)b_total;
            assert(size <= 1 && size >= 0);
            assert(m <= m_total && 0 <= m && m_total > 0);
            return 1 - fabs(m / m_total - size);
        };
    }

    auto linear_f(double a, double b) -> std::function<double(double, double, double, double)> {
        /*
         * Useful for finding a region of a certain size.
         */

        return [a, b] (double mv, double m_total, double bv, double b_total) {
            return a * mv / m_total + b * bv / b_total;
        };
    }

    std::function<double(double, double)> rho_f(std::function<double(double, double, double)> const& f, double rho) {
       return [&](double x, double y) {
           return f(x, y, rho);
       };
    }

    Subgrid maxSubgridLin(Grid const& grid, double eps, discrepancy_func_t const& f) {
      return max_subgrid_convex(grid, eps, f);
    }


    Subgrid maxSubgridSlow(Grid const &grid, discrepancy_func_t const& f) {
        return max_subgrid(grid, f);
    }


    double evaluate(discrepancy_func_t const& f, double m, double m_tot, double b, double b_tot) {
        return f(m, m_tot, b, b_tot);
    }

};


double evaluate_halfplane(pyscan::halfspace2_t const& d1, pyscan::wpoint_list_t const&  mpts, pyscan::wpoint_list_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_disk(pyscan::Disk const& d1, pyscan::wpoint_list_t const&  mpts, pyscan::wpoint_list_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_rectangle(pyscan::Rectangle const& d1, pyscan::wpoint_list_t const&  mpts, pyscan::wpoint_list_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_disk_labeled(pyscan::Disk const& d1, pyscan::lpoint_list_t const&  mpts, pyscan::lpoint_list_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_rectangle_labeled(pyscan::Rectangle const& d1, pyscan::lpoint_list_t const&  mpts, pyscan::lpoint_list_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_disk_traj(pyscan::Disk const& d1, pyscan::trajectory_set_t const&  mpts, pyscan::trajectory_set_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_rectangle_traj(pyscan::Rectangle const& d1, pyscan::trajectory_set_t const&  mpts, pyscan::trajectory_set_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_halfplane_traj(pyscan::halfspace2_t const& d1, pyscan::trajectory_set_t const&  mpts, pyscan::trajectory_set_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}

double evaluate_halfplane_labeled(pyscan::halfspace2_t const& d1, pyscan::lpoint_list_t const&  mpts, pyscan::lpoint_list_t const& bpts, pyscan::discrepancy_func_t const& disc) {
    return pyscan::evaluate_range(d1, mpts, bpts, disc);
}



PYBIND11_MODULE(libpyscan, pyscan_module){
    namespace py = pybind11;



    pyscan_module.def("hull", pyscan::graham_march);

    py::class_<pyscan::Grid>(pyscan_module, "Grid")
            .def(py::init<size_t, pyscan::wpoint_list_t const&, pyscan::wpoint_list_t const&>())
            .def("totalRedWeight", &pyscan::Grid::totalRedWeight)
            .def("totalBlueWeight", &pyscan::Grid::totalBlueWeight)
            .def("redWeight", &pyscan::Grid::redWeight)
            .def("blueWeight", &pyscan::Grid::blueWeight)
            .def("redSubWeight", &pyscan::Grid::redSubWeight)
            .def("blueSubWeight", &pyscan::Grid::blueSubWeight)
            .def("xCoord", &pyscan::Grid::xCoord)
            .def("yCoord", &pyscan::Grid::yCoord)
            .def("toRectangle", &pyscan::Grid::toRectangle)
            .def("size", &pyscan::Grid::size);

    py::class_<pyscan::Rectangle>(pyscan_module, "Rectangle")
            .def(py::init<double, double, double, double>())
            .def("lowX", &pyscan::Rectangle::lowX)
            .def("upX", &pyscan::Rectangle::upX)
            .def("lowY", &pyscan::Rectangle::lowY)
            .def("upY", &pyscan::Rectangle::upY)
            .def("__str__", &pyscan::Rectangle::toString)
            .def("contains", &pyscan::Rectangle::contains)
            .def("intersects_segment", &pyscan::Rectangle::intersects_segment)
            .def("intersects_trajectory", &pyscan::Rectangle::intersects_trajectory<bool>);

    py::class_<pyscan::Subgrid>(pyscan_module, "Subgrid")
            .def(py::init<size_t, size_t, size_t, size_t, double>())
            .def("lowCol", &pyscan::Subgrid::lowX)
            .def("upCol", &pyscan::Subgrid::upX)
            .def("lowRow", &pyscan::Subgrid::lowY)
            .def("upRow", &pyscan::Subgrid::upY)
            .def("__str__", &pyscan::Subgrid::toString)
            .def("__repr__", &pyscan::Subgrid::toString)
            .def("fValue", &pyscan::Subgrid::fValue);

    py::class_<pyscan::Disk>(pyscan_module, "Disk")
            .def(py::init<double, double, double>())
            .def("get_origin", &pyscan::Disk::getOrigin)
            .def("get_radius", &pyscan::Disk::getRadius)
            .def("contains", &pyscan::Disk::contains)
            .def("intersects_segment", &pyscan::Disk::intersects_segment)
            .def("__str__", &pyscan::Disk::str)
            .def("__repr__", &pyscan::Disk::str)
            .def("intersects_trajectory", &pyscan::Disk::intersects_trajectory<bool>);

    py::class_<pyscan::HalfSpace<2>>(pyscan_module, "Halfplane")
            .def(py::init<pyscan::Point<2>>())
            .def("get_coords", &pyscan::HalfSpace<2>::get_coords)
            .def("contains", &pyscan::HalfSpace<2>::contains)
            .def("intersects_segment", &pyscan::HalfSpace<2>::intersects_segment)
            .def("__str__", &pyscan::HalfSpace<2>::str)
            .def("__repr__", &pyscan::HalfSpace<2>::str)
            .def("intersects_trajectory", &pyscan::HalfSpace<2>::intersects_trajectory<bool>);

    py::class_<pyscan::HalfSpace<3>>(pyscan_module, "Halfspace")
            .def(py::init<pyscan::Point<3>>())
            .def("get_coords", &pyscan::HalfSpace<3>::get_coords)
            .def("contains", &pyscan::HalfSpace<3>::contains)
            .def("__str__", &pyscan::HalfSpace<3>::str)
            .def("__repr__", &pyscan::HalfSpace<3>::str)
            .def("intersects_segment", &pyscan::HalfSpace<3>::intersects_segment);

    py::class_<pyscan::pt2_t>(pyscan_module, "Point")
            .def(py::init<double, double, double>())
            .def(py::init<std::tuple<double, double>>())
            .def("approx_eq", &pyscan::Point<2>::approx_eq)
            .def("__getitem__", &pyscan::Point<2>::operator())
            .def("get_coord", &pyscan::Point<2>::get_coord)
            .def("above", &pyscan::Point<2>::above)
            .def("above_closed", &pyscan::Point<2>::above_closed)
            .def("below_closed", &pyscan::Point<2>::below_closed)
            .def("below", &pyscan::Point<2>::below)
            .def("crosses", &pyscan::Point<2>::crosses)
            .def("evaluate", &pyscan::Point<2>::evaluate)
            .def("orient_down", &pyscan::Point<2>::orient_down)
            .def("orient_up", &pyscan::Point<2>::orient_up)
            .def("dist", &pyscan::Point<2>::dist)
            .def("pdot", &pyscan::Point<2>::pdot)
            .def("parallel_lt", &pyscan::Point<2>::parallel_lt)
            .def("parallel_lte", &pyscan::Point<2>::parallel_lte)
            .def("__str__", &pyscan::Point<>::str)
            .def("__repr__", &pyscan::Point<>::str)
            .def("__eq__", &pyscan::Point<>::operator==);

    py::class_<pyscan::pt3_t>(pyscan_module, "Point3")
            .def(py::init<double, double, double, double>())
            .def(py::init<std::tuple<double, double, double>>())
            .def("approx_eq", &pyscan::Point<3>::approx_eq)
            .def("__getitem__", &pyscan::Point<3>::operator())
            .def("get_coord", &pyscan::Point<3>::get_coord)
            .def("above", &pyscan::Point<3>::above)
            .def("above_closed", &pyscan::Point<3>::above_closed)
            .def("below_closed", &pyscan::Point<3>::below_closed)
            .def("below", &pyscan::Point<3>::below)
            .def("crosses", &pyscan::Point<3>::crosses)
            .def("evaluate", &pyscan::Point<3>::evaluate)
            .def("orient_down", &pyscan::Point<3>::orient_down)
            .def("orient_up", &pyscan::Point<3>::orient_up)
            .def("dist", &pyscan::Point<3>::dist)
            .def("pdot", &pyscan::Point<3>::pdot)
            .def("parallel_lt", &pyscan::Point<3>::parallel_lt)
            .def("parallel_lte", &pyscan::Point<3>::parallel_lte)
            .def("__str__", &pyscan::Point<3>::str)
            .def("__repr__", &pyscan::Point<3>::str)
            .def("__eq__", &pyscan::Point<3>::operator==);

    /////////////////////////////////////////////////////////////////////
    //Kernel Scanning Wrappers//////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    py::class_<pyscan::kernel_func_t >(pyscan_module, "Kernel_f");


//    py::class_<pyscan::KDisc>(pyscan_module, "KDisc");
//
//    py::class_<pyscan::Bernouli_kf, pyscan::KDisc>(pyscan_module, "Bernoulli_K")
//            .def(py::init<pyscan::kernel_func_t, double>());
//
//
//    pyscan_module.def("max_annuli", &pyscan::max_annuli);
//    pyscan_module.def("max_annuli_scale", &pyscan::max_annuli_scale);
//    pyscan_module.def("max_annuli_scale_multi", &pyscan::max_annuli_scale_multi);

//    py::class_<pyscan::Bernoulli_Disk>(pyscan_module, "Bernoulli");
    pyscan_module.def("max_kernel", &pyscan::max_kernel);
    pyscan_module.def("max_kernel_prune_far", &pyscan::max_kernel_prune_far);
    pyscan_module.def("max_kernel_adaptive", &pyscan::max_kernel_adaptive);
    pyscan_module.def("max_kernel_slow2", &pyscan::max_kernel_slow2);

    pyscan_module.def("max_kernel_slow", &pyscan::max_kernel_slow);
    pyscan_module.def("measure_kernel", &pyscan::measure_kernel);

//    pyscan_module.attr("GAUSSIAN_KERNEL") = pyscan::kernel_func_t(
//            [](double dist, double bandwidth) {
//                return exp(- dist * dist / (bandwidth * bandwidth));
//    });


    ////////////////////////////////////////////////////////////////////
    //TrajectoryScan.hpp wrappers///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    py::class_<pyscan::trajectory_t>(pyscan_module, "Trajectory")
            .def(py::init<pyscan::point_list_t>())
            .def(py::init<std::vector<std::tuple<double, double>>>())
            .def("point_dist", &pyscan::trajectory_t::point_dist)
            .def("get_weight", &pyscan::trajectory_t::get_weight)
            .def("get_length", &pyscan::trajectory_t::get_length)
            .def("point_dist", &pyscan::trajectory_t::point_dist)
            .def("get_pts", &pyscan::trajectory_t::get_pts);

    py::class_<pyscan::wtrajectory_t>(pyscan_module, "WTrajectory")
            .def(py::init<double, pyscan::point_list_t>())
            .def(py::init<double, std::vector<std::tuple<double, double>>>())
            .def("point_dist", &pyscan::wtrajectory_t::point_dist)
            .def("get_weight", &pyscan::wtrajectory_t::get_weight)
            .def("get_length", &pyscan::wtrajectory_t::get_length)
            .def("point_dist", &pyscan::wtrajectory_t::point_dist)
            .def("get_pts", &pyscan::wtrajectory_t::get_pts);


    py::implicitly_convertible<std::tuple<double, double> , pyscan::Point<2> >();
    py::implicitly_convertible<std::tuple<double, double, double> , pyscan::Point<3> >();
    py::implicitly_convertible<pyscan::point_list_t, pyscan::Trajectory>();

    
    pyscan_module.def("to_segment", &pyscan::to_Segment);

    pyscan_module.def("aeq", &pyscan::util::aeq);

    py::class_<pyscan::Segment, pyscan::HalfSpace<2> >(pyscan_module, "Segment")
            .def(py::init<pyscan::HalfSpace<2>, pyscan::Point<>, pyscan::Point<>>())
            .def("lte", &pyscan::Segment::lte)
            .def("lt", &pyscan::Segment::lt)
            .def("gt", &pyscan::Segment::gt)
            .def("gte", &pyscan::Segment::gte)
            .def("crossed", &pyscan::Segment::crossed)
            .def("split", &pyscan::Segment::split)
//            .def("__str__", &pyscan::Segment::str)
//            .def("__repr__", &pyscan::Segment::str)
            .def("get_e1", &pyscan::Segment::get_e1)
            .def("get_e2", &pyscan::Segment::get_e2);

    py::class_<pyscan::WPoint<2>, pyscan::Point<2>>(pyscan_module, "WPoint")
        .def(py::init<double, double, double, double>())
        .def("get_weight", &pyscan::WPoint<2>::get_weight);

    py::class_<pyscan::LPoint<2>, pyscan::WPoint<2>>(pyscan_module, "LPoint")
        .def(py::init<size_t, double, double, double, double>())
        .def("get_label", &pyscan::LPoint<2>::get_label);

    py::class_<pyscan::WPoint<3>, pyscan::Point<3>>(pyscan_module, "WPoint3")
        .def(py::init<double, double, double, double, double>())
        .def("get_weight", &pyscan::WPoint<3>::get_weight);

    py::class_<pyscan::LPoint<3>, pyscan::WPoint<3>>(pyscan_module, "LPoint3")
        .def(py::init<size_t, double, double, double, double, double>())
        .def("get_label", &pyscan::LPoint<3>::get_label);

    py::class_<pyscan::discrepancy_func_t >(pyscan_module, "CFunction");

    pyscan_module.attr("KULLDORF") = pyscan::discrepancy_func_t(
        [&](double m, double m_tot, double b, double b_tot) {
            return pyscan::kulldorff(m / m_tot, b / b_tot, .0001);
    });

    pyscan_module.attr("DISC") = pyscan::discrepancy_func_t (
        [&](double m, double m_tot, double b, double b_tot) {
            return std::abs(m / m_tot - b / b_tot);
    });

    pyscan_module.attr("RKULLDORF") = pyscan::discrepancy_func_t(
        [&](double m, double m_tot, double b, double b_tot) {
            return pyscan::regularized_kulldorff(m / m_tot, b / b_tot, .0001);
    });

    pyscan_module.def("bernoulli", pyscan::get_bernoulli);
    pyscan_module.def("rbernoulli", pyscan::get_bernoulli_single_sample);

    pyscan_module.def("evaluate", &pyscan::evaluate);
    pyscan_module.def("size_region", &pyscan::sized_region);
    pyscan_module.def("linear_f", &pyscan::linear_f);

    pyscan_module.def("intersection", &pyscan::intersection);
    pyscan_module.def("correct_orientation", &pyscan::correct_orientation);


    pyscan_module.def("max_subgrid", &pyscan::max_subgrid);
    pyscan_module.def("max_subgrid_convex", &pyscan::max_subgrid_convex);
    pyscan_module.def("max_subgrid_linear", &pyscan::max_subgrid_linear);
    pyscan_module.def("max_rectangle", &pyscan::max_rectangle);

    pyscan_module.def("make_net_grid", &pyscan::make_net_grid);
    pyscan_module.def("make_exact_grid", &pyscan::make_exact_grid);

    //Max Halfspace codes
    pyscan_module.def("max_halfplane", &pyscan::max_halfplane);
    pyscan_module.def("max_halfplane_labeled", &pyscan::max_halfplane_labeled);
    pyscan_module.def("max_halfspace", &pyscan::max_halfspace);
    pyscan_module.def("max_halfspace_labeled", &pyscan::max_halfspace_labeled);
    //pyscan_module.def("max_halfplane_fast", &pyscan::max_halfplane_fast);
    pyscan_module.def("ham_tree_sample", &pyscan::ham_tree_sample);

    pyscan_module.def("max_disk", &pyscan::max_disk);
    pyscan_module.def("max_disk_labeled", &pyscan::max_disk_labeled);
    //pyscan_module.def("max_disk_lift_labeled", &pyscan::max_disk_labeled);


    pyscan_module.def("evaluate_disk", &evaluate_disk);
    pyscan_module.def("evaluate_disk_alt", &evaluate_disk);
    pyscan_module.def("evaluate_disk_labeled", &evaluate_disk_labeled);
    pyscan_module.def("evaluate_disk_trajectory", &evaluate_disk_traj);

    pyscan_module.def("evaluate_halfplane", &evaluate_halfplane);
    pyscan_module.def("evaluate_halfplane_labeled", &evaluate_halfplane_labeled);
    pyscan_module.def("evaluate_halfplane_trajectory", &evaluate_halfplane_traj);

    pyscan_module.def("evaluate_rectangle", &evaluate_rectangle);
    pyscan_module.def("evaluate_rectangle_labeled", &evaluate_rectangle_labeled);
    pyscan_module.def("evaluate_rectangle_trajectory", &evaluate_rectangle_traj);


//    pyscan_module.def("max_disk_cached", &pyscan::max_disk_cached);
//    pyscan_module.def("max_disk_cached_labeled", &pyscan::max_disk_cached_labeled);

    pyscan_module.def("max_disk_scale", &pyscan::max_disk_scale);
    pyscan_module.def("max_disk_scale_labeled", &pyscan::max_disk_scale_labeled);
    pyscan_module.def("max_disk_scale_labeled_alt", &pyscan::max_disk_scale_labeled_alt);

    pyscan_module.def("max_rect_labeled", &pyscan::max_rect_labeled);
    pyscan_module.def("max_rectangle", &pyscan::max_rectangle);


    pyscan_module.def("max_rect_labeled_scale", pyscan::max_rect_labeled_scale);


    pyscan_module.def("max_disk_traj_grid", &pyscan::max_disk_traj_grid);
//

    //This simplifies the trajectory by using the dp algorithm.
    pyscan_module.def("dp_compress", &pyscan::dp_compress);
    //This grids the trajectory and assigns a single point to each cell.
    pyscan_module.def("grid_kernel", &pyscan::approx_traj_grid);
    pyscan_module.def("grid_trajectory", &pyscan::grid_traj);
    //This grids the trajectory and creates an alpha hull in each one.
    pyscan_module.def("grid_direc_kernel", &pyscan::approx_traj_kernel_grid);
    //This is for 2d eps-kernel useful for halfspaces.
    pyscan_module.def("halfplane_kernel", pyscan::approx_hull);
    pyscan_module.def("convex_hull", pyscan::graham_march);
    //This is a 3d eps-kernel for disks.
    pyscan_module.def("lifting_kernel", &pyscan::lifting_coreset);

    pyscan_module.def("coreset_error_halfplane", &pyscan::error_halfplane_coreset);
    //pyscan_module.def("coreset_error_disk", &pyscan::error_disk_coreset);

    //This is for partial scanning, but could be used for full scannings.
    pyscan_module.def("block_sample", &pyscan::block_sample);
    pyscan_module.def("uniform_sample", &pyscan::uniform_sample);
    pyscan_module.def("even_sample", &pyscan::even_sample);

    pyscan_module.def("block_sample_error", &pyscan::block_sample_error);
    pyscan_module.def("uniform_sample_error", &pyscan::uniform_sample_error);
    pyscan_module.def("even_sample_error", &pyscan::even_sample_error);

    pyscan_module.def("polygon_sample", &pyscan::polygon_sample);
    pyscan_module.def("polygon_grid", &pyscan::polygon_grid);

    pyscan_module.def("polygon_grid_even", &pyscan::polygon_grid_even);
    pyscan_module.def("polygon_grid_hull", &pyscan::polygon_grid_hull);

    pyscan_module.def("naive_find_rect", &pyscan::naive_find_rect);
    pyscan_module.def("naive_scan_grid", &pyscan::naive_scan_grid);
    pyscan_module.def("naive_approx_find_rect", &pyscan::naive_approx_find_rect);
    pyscan_module.def("scan_grid", &pyscan::scan_grid);
    pyscan_module.def("find_rect", &pyscan::find_rect);


    //Satscan comparison function
    pyscan_module.def("satscan_grid", pyscan::satscan_grid);
    pyscan_module.def("satscan_labeled", pyscan::satscan_grid_labeled);

    pyscan_module.def("kernel_centers", pyscan::kernel_centers_approximate);

}
