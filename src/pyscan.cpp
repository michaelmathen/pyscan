//
// Created by mmath on 7/7/17.
//

#include "RectangleScan.hpp"
#include "EpsSamples.hpp"
#include "HalfplaneScan.hpp"
#include "DiskScan.hpp"

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

auto toLine(py::object const& el)-> pyscan::Line {
    if (py::len(el) == 2) {
        return pyscan::Line(py::extract<double>(el[0]), py::extract<double>(el[1]), 0);
    } else {
        return pyscan::Line(py::extract<double>(el[0]), py::extract<double>(el[1]), py::extract<int>(el[2]));
    }
}

auto toPoint(py::object const& el)-> pyscan::Line2 {
    return pyscan::Line2(py::extract<double>(el[0]), py::extract<double>(el[1]));
}

auto toLineList(py::object const& el) -> pyscan::LineList {
    pyscan::LineList l;
    for (auto beg = py::stl_input_iterator<py::object>(el);
         beg != py::stl_input_iterator<py::object>(); beg++) {
        l.push_back(toLine(*beg));
    }
    return l;
}

auto toPointList(py::object const& el) -> pyscan::PointList {
    pyscan::PointList l;
    for (auto beg = py::stl_input_iterator<py::object>(el);
         beg != py::stl_input_iterator<py::object>(); beg++) {
        l.push_back(toPoint(*beg));
    }
    return l;
}

py::tuple toPyLine2(pyscan::Line const& line) {
    return py::make_tuple(std::get<0>(line), std::get<1>(line));
}

py::tuple toPyPoint(pyscan::Line2 const& pt) {
    return py::make_tuple(std::get<0>(pt), std::get<1>(pt));
}

template <typename T, typename F>
py::list toPyList(T const& lists, F converter) {
    py::list l;
    for (auto& point : lists) {
        l.append(converter(point));
    }
    return l;
}

py::list cellsToList(const std::vector<pyscan::Cell>& v) {
    py::list l;
    for (auto cell : v) {
        l.append(py::make_tuple(std::get<0>(cell),
                                toPyList(std::get<1>(cell), toPyLine2),
                                toPyList(std::get<2>(cell), toPyPoint),
                                std::get<3>(cell)));
    }
    return l;
}

py::list testSetWrapper(py::object& iterable, int r) {
    pyscan::LineList lines;
    py::stl_input_iterator<py::object> end;
    for (auto obj = py::stl_input_iterator<py::object>(iterable); end != obj; obj++) {
        lines.push_back(toLine(*obj));
    }
    pyscan::PointList points;
    auto cells = pyscan::randomCuttings(lines, points, r);
    auto test_set = pyscan::dualizeToLines(cells);

    py::list output_lines;
    for (auto line : test_set) {
        output_lines.append(py::make_tuple(std::get<0>(line), std::get<1>(line)));
    }
    return output_lines;
}

py::list splitCellWrapper(pyscan::Trapezoid const& trap, py::object const& lines, py::object const& line) {
    pyscan::LineList lineVector = toLineList(lines);
    pyscan::PointList pointVector;
    std::vector<pyscan::Cell> output;
    pyscan::splitCell(pyscan::Cell(trap, lineVector, pointVector, lineVector.size()), toLine(line), output);
    return cellsToList(output);
}

pyscan::Grid makeGrid(const py::object& iterable, size_t r) {
    std::vector<pyscan::Point<>> points = to_std_vector<pyscan::Point<>>(iterable);

    return pyscan::Grid(points.begin(), points.end(), r);
}

pyscan::Grid makeNetGrid(const py::object& iterable, size_t r) {
    std::vector<pyscan::Point<>> points = to_std_vector<pyscan::Point<>>(iterable);

    return pyscan::Grid(points.begin(), points.end(), r, true);
}

pyscan::Halfspace<2> maxHalfplaneStat(const py::object& net, const py::object& sample, double rho) {
    std::vector<pyscan::Point<>> net_points = to_std_vector<pyscan::Point<>>(net);
    std::vector<pyscan::Point<>> sample_points = to_std_vector<pyscan::Point<>>(sample);
    auto pairs = maxHalfplaneStat(net_points.begin(), net_points.end(), sample_points.begin(), sample_points.end(), rho);
    return std::get<0>(pairs);
}

pyscan::Halfspace<2> maxHalfplaneGamma(const py::object& net, const py::object& sample, double rho) {
    std::vector<pyscan::Point<>> net_points = to_std_vector<pyscan::Point<>>(net);
    std::vector<pyscan::Point<>> sample_points = to_std_vector<pyscan::Point<>>(sample);
    auto pairs = maxHalfplaneGamma(net_points.begin(), net_points.end(), sample_points.begin(), sample_points.end(), rho);
    return std::get<0>(pairs);
}

pyscan::Halfspace<2> maxHalfplaneLin(const py::object& net, const py::object& sample, double rho) {
    std::vector<pyscan::Point<>> net_points = to_std_vector<pyscan::Point<>>(net);
    std::vector<pyscan::Point<>> sample_points = to_std_vector<pyscan::Point<>>(sample);
    auto pairs = maxHalfplaneLin(net_points.begin(), net_points.end(), sample_points.begin(), sample_points.end());
    return std::get<0>(pairs);
}


pyscan::Disk maxDisk(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    std::vector<pyscan::Point<>> net_points = to_std_vector<pyscan::Point<>>(net);
    std::vector<pyscan::Point<>> sample_p_M = to_std_vector<pyscan::Point<>>(sampleM);
    std::vector<pyscan::Point<>> sample_p_B = to_std_vector<pyscan::Point<>>(sampleB);

    return diskScanStat(net_points, sample_p_M, sample_p_B, rho);
}



pyscan::Disk maxDiskLabelsD(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    auto net_points = to_std_vector<pyscan::LPoint<>>(net);
    auto sample_p_M = to_std_vector<pyscan::LPoint<>>(sampleM);
    auto sample_p_B = to_std_vector<pyscan::LPoint<>>(sampleB);
    return diskScanStatLabels(net_points, sample_p_M, sample_p_B, rho);
}

/*
pyscan::Disk maxDiskLabelsI(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    auto net_points = to_std_vector<pyscan::LPoint<int>>(net);
    auto sample_p_M = to_std_vector<pyscan::LPoint<int>>(sampleM);
    auto sample_p_B = to_std_vector<pyscan::LPoint<int>>(sampleB);
    return diskScanStatLabels(net_points, sample_p_M, sample_p_B, rho);
}
 */
pyscan::Rectangle maxRectLabelsD(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    auto net_points = to_std_vector<pyscan::LPoint<>>(net);
    auto sample_p_M = to_std_vector<pyscan::LPoint<>>(sampleM);
    auto sample_p_B = to_std_vector<pyscan::LPoint<>>(sampleB);
    return pyscan::maxRectStatLabels(net_points, sample_p_M, sample_p_B, rho);
}
/*
pyscan::Subgrid maxRectLabelsI(const py::object& net, const py::object& sampleM, const py::object& sampleB, double rho) {
    auto net_points = to_std_vector<pyscan::LPoint<int>>(net);
    auto sample_p_M = to_std_vector<pyscan::LPoint<int>>(sampleM);
    auto sample_p_B = to_std_vector<pyscan::LPoint<int>>(sampleB);
    return pyscan::maxLabeledRectStat(net_points, sample_p_M, sample_p_B, rho);
}
*/
double local_xLoc(double a1, double a2, double a3, double a4) {
    return pyscan::xLoc(a1, a2, a3, a4);
}

py::list randomCuttings(const py::object& iterable, int r) {
    pyscan::LineList lines;
    py::stl_input_iterator<py::object> end;
    for (auto obj = py::stl_input_iterator<py::object>(iterable); end != obj; obj++) {
        lines.push_back(toLine(*obj));
    }
    pyscan::PointList points;
    return cellsToList(pyscan::randomCuttings(lines, points, r));
}

py::list createPartitionWrapper(const py::object& iterable, int r) {
    pyscan::PointList points;
    py::stl_input_iterator<py::object> end;
    for (auto obj = py::stl_input_iterator<py::object>(iterable); end != obj; obj++) {
        points.push_back(toPoint(*obj));
    }
    return cellsToList(pyscan::createPartition(points, r));
}

BOOST_PYTHON_MODULE(pyscan) {
    using namespace py;

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

    py::class_<pyscan::Point<2>>("Point", py::init<double, double, double, double>())
            .def("getWeight", &pyscan::Point<>::getWeight)
            .def("setRedWeight", &pyscan::Point<>::setRedWeight)
            .def("setBlueWeight", &pyscan::Point<>::setBlueWeight)
            .def("getRedWeight", &pyscan::Point<>::getRedWeight)
            .def("getBlueWeight", &pyscan::Point<>::getBlueWeight)
            .def("__str__", &pyscan::Point<>::toString)
            .def("__repr__", &pyscan::Point<>::toString)
            .def("getX", &pyscan::Point<>::getX)
            .def("getY", &pyscan::Point<>::getY);

    /*
    py::class_<pyscan::Point<double, 3>>("__point3d", py::init<double, double, double, double, double>())
            .def("getWeight", &pyscan::Point<double, 3>::getWeight)
	          .def("setRedWeight", &pyscan::Point<double, 3>::setRedWeight)
	          .def("setBlueWeight", &pyscan::Point<double, 3>::setBlueWeight)
            .def("getRedWeight", &pyscan::Point<double, 3>::getRedWeight)
            .def("getBlueWeight", &pyscan::Point<double, 3>::getBlueWeight)
            .def("__str__", &pyscan::Point<double, 3>::toString)
            .def("__repr__", &pyscan::Point<double, 3>::toString);
    */
    py::class_<pyscan::Disk>("Disk", py::init<double, double, double, double>())
            .add_property("a", &pyscan::Disk::getA, &pyscan::Disk::setA)
            .add_property("b", &pyscan::Disk::getB, &pyscan::Disk::setB)
            .add_property("radius", &pyscan::Disk::getR, &pyscan::Disk::setR)
            .def("contains", &pyscan::Disk::contains)
            .def("fValue", &pyscan::Disk::fValue);

    py::class_<pyscan::LPoint<2>, py::bases<pyscan::Point<2>>>("LPoint", py::init<size_t, double, double, double, double>())
	          .def("getLabel", &pyscan::LPoint<2>::getLabel);

    py::class_<pyscan::BloomFilter>("BloomFilter", py::init<int, double>())
            .def("insert", &pyscan::BloomFilter::insert)
            .def("mightBePresent", &pyscan::BloomFilter::mightBePresent);

    py::class_<pyscan::Trapezoid>("Trapezoid", py::init<double, double, double, double, double, double>())
            .def("crossing", &pyscan::Trapezoid::crossing)
            .def("crossesLeftSide", &pyscan::Trapezoid::crossesLeftSide)
            .def("crossesRightSide", &pyscan::Trapezoid::crossesRightSide)
            .def("crossesTop", &pyscan::Trapezoid::crossesTop)
            .def("crossesBottom", &pyscan::Trapezoid::crossesBottom)
            .def_readonly("top_a", &pyscan::Trapezoid::top_a)
            .def_readonly("top_b", &pyscan::Trapezoid::top_b)
            .def_readonly("left_x", &pyscan::Trapezoid::left_x)
            .def_readonly("right_x", &pyscan::Trapezoid::right_x)
            .def_readonly("bottom_a", &pyscan::Trapezoid::bottom_a)
            .def("__print__", &pyscan::Trapezoid::print)
            .def("__repr__", &pyscan::Trapezoid::print)
            .def_readonly("bottom_b", &pyscan::Trapezoid::bottom_b);

    py::class_<pyscan::Halfspace<2>>("__halfspace2", py::init<double, double, double>())
            .def("geta", &pyscan::Halfspace<2>::get<0>)
            .def("getb", &pyscan::Halfspace<2>::get<1>)
            .def("fValue", &pyscan::Halfspace<2>::fValue);

    py::class_<pyscan::Halfspace<3>>("__halfspace3", py::init<double, double, double, double>())
            .def("geta1", &pyscan::Halfspace<3>::get<0>)
            .def("geta2", &pyscan::Halfspace<3>::get<1>)
            .def("getb", &pyscan::Halfspace<3>::get<2>)
            .def("fValue", &pyscan::Halfspace<3>::fValue);


    /*
     class_<pyscan::SlabTree<int>>("SlabTree", py::init<pyscan::Grid<int>, int>())
            .def("measure", &pyscan::SlabTree<int>::measure)
            .def("measureSubgrid", &pyscan::SlabTree<int>::measureSubgrid);
    */

    py::def("testSet", testSetWrapper);
    py::def("createPartition", createPartitionWrapper);
    py::def("xLoc", local_xLoc);
    py::def("splitCell", splitCellWrapper);
    py::def("randomCuttings", randomCuttings);
    py::def("makeGrid", makeGrid);
    py::def("makeNetGrid", makeNetGrid);


    pyscan::Subgrid (*f1)(pyscan::Grid const&, double) = &pyscan::maxSubgridKullSlow;
    pyscan::Subgrid (*f2)(pyscan::Grid const&, double, double) = &pyscan::maxSubgridLinearSlow;

    py::def("maxSubgridKullSlow", f1);
    py::def("maxSubgridLinearSlow", f2);

    pyscan::Subgrid (*m1)(pyscan::Grid const&, double, double) = &pyscan::maxSubgridLinearSimple;
    py::def("maxSubgridLinearSimple", m1);
    py::def("maxSubgridLinear", &pyscan::maxSubgridLinearG);

    pyscan::Subgrid (*fl1)(pyscan::Grid const&, double) = &pyscan::maxSubgridLinKull;
    py::def("maxSubgridKull", fl1);

    py::def("maxHalfPlaneStat", &maxHalfplaneStat);
    py::def("maxHalfPlaneLin", &maxHalfplaneLin);
    py::def("maxHalfPlaneGamma", &maxHalfplaneGamma);
    py::def("maxDisk", &maxDisk);
    //py::def("maxDiskLabels", &maxDiskLabelsI);
    py::def("maxDiskLabels", &maxDiskLabelsD);
    //py::def("maxRectLabels", &maxRectLabelsI);
    py::def("maxRectLabels", &maxRectLabelsD);


}
