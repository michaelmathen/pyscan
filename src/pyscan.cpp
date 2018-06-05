//
// Created by mmath on 7/7/17.
//

#include <functional>

#include "RectangleScan.hpp"
#include "HalfplaneScan.hpp"
#include "DiskScan.hpp"
#include "DiskScan2.hpp"

#include <boost/iterator_adaptors.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <iostream>

namespace py = boost::python;


template<typename T>
std::vector<T> to_std_vector(const py::object& iterable) {
    return std::vector<T>(py::stl_input_iterator<T>(iterable),
                          py::stl_input_iterator<T>());
}

auto toPointList(py::object const& el) -> pyscan::point_list {
    pyscan::point_list l;
    for (auto beg = py::stl_input_iterator<pyscan::Point<>>(el);
         beg != py::stl_input_iterator<pyscan::Point<>>(); beg++) {
        l.push_back(*beg);
    }
    return l;
}




namespace pyscan {



    auto toLPt(const py::object& pts, 
            const py::object& weights,
            const py::object& labels) -> pyscan::lpoint_list {

        auto bp = py::stl_input_iterator<Pt2>(pts);
        auto ep = py::stl_input_iterator<Pt2>();
        auto bw = py::stl_input_iterator<double>(pts);
        auto ew = py::stl_input_iterator<double>();
        auto bl = py::stl_input_iterator<size_t>(labels);
        auto el = py::stl_input_iterator<size_t>();
        //assert(((ep - bp) == (ew - bw)) && ((ew -bw) == (el - bl)));

        pyscan::lpoint_list list;
        while(bp != ep) {
            list.emplace_back(pyscan::LPoint<>(*bl, *bw, (*bp)[0], (*bp)[1], (*bp)[2]));
            bp++;
            bw++;
            bl++;
        }
        return list;
    }

    auto toWPt(const py::object& pts, 
            const py::object& weights) -> pyscan::wpoint_list {

        auto bp = py::stl_input_iterator<Pt2>(pts);
        auto ep = py::stl_input_iterator<Pt2>();
        auto bw = py::stl_input_iterator<double>(pts);
        auto ew = py::stl_input_iterator<double>();
        //assert(((ep - bp) == (ew - bw)) && ((ew -bw) == (el - bl)));

        pyscan::wpoint_list list;
        while(bp != ep) {
            list.emplace_back(pyscan::WPoint<>(*bw, (*bp)[0], (*bp)[1], (*bp)[2]));
            bp++;
            bw++;
        }
        return list;
    }

    std::function<double(double, double)> rho_f(std::function<double(double, double, double)> const& f, double rho) {
       return [&](double x, double y) {
           return f(x, y, rho);
       };
    }

    Subgrid maxSubgridLinKull(Grid const& grid, double eps, double rho) {
      return maxSubgridLinearSimple(grid, eps, rho_f(regularized_kulldorff, rho));
    }

    Subgrid maxSubgridLinKullTheory(Grid const& grid, double eps, double rho){
        return maxSubgridLinearTheory(grid, eps, rho_f(regularized_kulldorff, rho));
    }

    Subgrid maxSubgridLinGamma(Grid const& grid, double eps) {
        return maxSubgridLinearSimple(grid, eps, rho_f(gamma, 0.0));
    }

    Subgrid maxSubgridKullSlow(Grid const &grid, double rho) {
        return maxSubgridNonLinear(grid, rho_f(kulldorff, rho));
    }


    Subgrid maxSubgridLinearSlow(Grid const& grid, double a, double b) {
        return maxSubgridNonLinear(grid, [&](double red, double blue) {
            return a * red + b * blue;
        });
    }

    py::tuple maxDisk(const py::object& net, 
            const py::object& sampleM, 
            const py::object& weightM,
            const py::object& sampleB,
            const py::object& weightB,
             double rho) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto sample_p_M = toWPt(sampleM, weightM);
        auto sample_p_B = toWPt(sampleB, weightB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = diskScanSlow(net_points, sample_p_M, 
            sample_p_B, 
            rho_f(kulldorff, rho));
        return py::make_tuple(d1, d1value);
    }

    py::tuple maxDiskLabels(const py::object& net, 
            const py::object& sampleM, 
            const py::object& weightM,
            const py::object& labelM,
            const py::object& sampleB,
            const py::object& weightB,
            const py::object& labelB,
             double rho) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto lpointsM = toLPt(sampleM, weightM, labelM);
        auto lpointsB = toLPt(sampleB, weightB, labelB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = diskScanLabels(net_points, lpointsM, lpointsB, rho_f(kulldorff, rho));
        return py::make_tuple(d1, d1value);
    }

    // pyscan::Rectangle maxRectLabelsD(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    //     auto net_points = to_std_vector<pyscan::LPoint<>>(net);
    //     auto sample_p_M = to_std_vector<pyscan::LPoint<>>(sampleM);
    //     auto sample_p_B = to_std_vector<pyscan::LPoint<>>(sampleB);
    //     return pyscan::maxRectStatLabels(net_points, sample_p_M, sample_p_B, rho_f(kulldorff, rho));
    // }
};





pyscan::Grid makeGrid(const py::object& sample_r, const py::object& weight_r, const py::object& sample_b,
                      const py::object& weight_b, size_t r) {
    auto points_r = to_std_vector<pyscan::Point<>>(sample_r);
    auto points_b = to_std_vector<pyscan::Point<>>(sample_b);
    auto v_weight_r = to_std_vector<double>(weight_r);
    auto v_weight_b = to_std_vector<double>(weight_b);
    return pyscan::Grid(r, points_r, v_weight_r, points_b, v_weight_b);
}

//pyscan::Grid makeNetGrid(const py::object& net, const py::object& sample_r, const py::object& sample_b) {
//    std::vector<pyscan::Point<>> net_p = to_std_vector<pyscan::Point<>>(net);
//    std::vector<pyscan::Point<>> points_r = to_std_vector<pyscan::Point<>>(sample_r);
//    std::vector<pyscan::Point<>> points_b = to_std_vector<pyscan::Point<>>(sample_b);
//    return pyscan::Grid(net_p.begin(), net_p.end(), points_r.begin(), points_r.end(), points_b.begin(), points_b.end());
//}


/*
pyscan::Disk maxDiskLabelsI(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    auto net_points = to_std_vector<pyscan::LPoint<int>>(net);
    auto sample_p_M = to_std_vector<pyscan::LPoint<int>>(sampleM);
    auto sample_p_B = to_std_vector<pyscan::LPoint<int>>(sampleB);
    return diskScanStatLabels(net_points, sample_p_M, sample_p_B, rho);
}
 */


BOOST_PYTHON_MODULE(pyscan) {
    using namespace py;
    /*
    py::class_<pyscan::MaximumIntervals>("MaximumIntervals", py::init<std::size_t>())
            .def("mergeZeros")
    */
    py::class_<pyscan::Grid>("Grid", py::no_init)
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

    py::class_<pyscan::Rectangle>("Rectangle", py::init<double, double, double, double, double>())
            .def("lowX", &pyscan::Rectangle::lowX)
            .def("upX", &pyscan::Rectangle::upX)
            .def("lowY", &pyscan::Rectangle::lowY)
            .def("upY", &pyscan::Rectangle::upY)
            .def("__str__", &pyscan::Rectangle::toString)
            .def("contains", &pyscan::Rectangle::contains)
            .def("fValue", &pyscan::Rectangle::fValue);

    py::class_<pyscan::Subgrid>("Subgrid", py::init<size_t, size_t, size_t, size_t, double>())
            .def("lowCol", &pyscan::Subgrid::lowX)
            .def("upCol", &pyscan::Subgrid::upX)
            .def("lowRow", &pyscan::Subgrid::lowY)
            .def("upRow", &pyscan::Subgrid::upY)
            .def("__str__", &pyscan::Subgrid::toString)
            .def("__repr__", &pyscan::Subgrid::toString)
            .def("fValue", &pyscan::Subgrid::fValue);

    py::class_<pyscan::Pt2>("Point", py::init<double, double>())
            .def("approx_eq", &pyscan::Point<>::approx_eq)
            .def("__getitem__", &pyscan::Point<>::operator[])
            .def("above", &pyscan::Point<>::above)
            .def("above_closed", &pyscan::Point<>::above_closed)
            .def("parallel", &pyscan::Point<>::parallel)
            .def("getX", &pyscan::getX<2>)
            .def("getY", &pyscan::getY<2>);

    py::def("dot", &pyscan::dot<2>);
    py::def("intersection", &pyscan::intersection);
    py::def("correct_orientation", &pyscan::correct_orientation);
    py::def("above_closed_interval", &pyscan::above_closed_interval);
    py::def("above_interval", &pyscan::above_interval);


    py::class_<pyscan::Disk>("Disk", py::init<double, double, double>())
            .add_property("a", &pyscan::Disk::getA, &pyscan::Disk::setA)
            .add_property("b", &pyscan::Disk::getB, &pyscan::Disk::setB)
            .add_property("radius", &pyscan::Disk::getR, &pyscan::Disk::setR)
            .def("contains", &pyscan::Disk::contains);


    py::def("makeGrid", makeGrid);
    //py::def("makeNetGrid", makeNetGrid);
    //py::def("makeSampleGrid", makeSampleGrid);

    pyscan::Subgrid (*f1)(pyscan::Grid const&, double) = &pyscan::maxSubgridKullSlow;
    pyscan::Subgrid (*f2)(pyscan::Grid const&, double, double) = &pyscan::maxSubgridLinearSlow;

    py::def("maxSubgridKullSlow", f1);
    py::def("maxSubgridLinearSlow", f2);

    pyscan::Subgrid (*m1)(pyscan::Grid const&, double, double) = &pyscan::maxSubgridLinearSimple;
    py::def("maxSubgridLinearSimple", m1);
    py::def("maxSubgridLinear", &pyscan::maxSubgridLinearG);

    pyscan::Subgrid (*fl1)(pyscan::Grid const&, double, double) = &pyscan::maxSubgridLinKull;
    py::def("maxSubgridKull", fl1);

    py::def("maxSubgridTheoryKull", pyscan::maxSubgridLinKullTheory);

    //py::def("maxDiskLabels", &maxDiskLabelsI);
    py::def("maxDiskLabels", &pyscan::maxDiskLabels);
    //py::def("maxRectLabels", &maxRectLabelsI);
    //py::def("maxRectLabels", &maxRectLabelsD);


}
