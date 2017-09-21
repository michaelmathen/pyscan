//
// Created by mmath on 7/7/17.
//

#include "RectangleScan.hpp"
#include "EpsSamples.hpp"
#include "HalfplaneScan.hpp"


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

pyscan::Grid<int> makeGrid(const py::object& iterable, size_t r) {
    std::vector<pyscan::Point<int>> points = to_std_vector<pyscan::Point<int>>(iterable);

    return pyscan::Grid<int>(points.begin(), points.end(), r);
}

pyscan::Grid<int> makeNetGrid(const py::object& iterable, size_t r) {
    std::vector<pyscan::Point<int>> points = to_std_vector<pyscan::Point<int>>(iterable);

    return pyscan::Grid<int>(points.begin(), points.end(), r, true);
}

pyscan::Halfplane maxHalfplaneStat(const py::object& net, const py::object& sample, double rho) {
    std::vector<pyscan::Point<double>> net_points = to_std_vector<pyscan::Point<double>>(net);
    std::vector<pyscan::Point<double>> sample_points = to_std_vector<pyscan::Point<double>>(sample);
    return maxHalfplaneStat(net_points.begin(), net_points.end(), sample_points.begin(), sample_points.end(), rho);
}

pyscan::Halfplane maxHalfplaneGamma(const py::object& net, const py::object& sample, double rho) {
    std::vector<pyscan::Point<double>> net_points = to_std_vector<pyscan::Point<double>>(net);
    std::vector<pyscan::Point<double>> sample_points = to_std_vector<pyscan::Point<double>>(sample);
    return maxHalfplaneGamma(net_points.begin(), net_points.end(), sample_points.begin(), sample_points.end(), rho);
}

pyscan::Halfplane maxHalfplaneLin(const py::object& net, const py::object& sample, double rho) {
    std::vector<pyscan::Point<double>> net_points = to_std_vector<pyscan::Point<double>>(net);
    std::vector<pyscan::Point<double>> sample_points = to_std_vector<pyscan::Point<double>>(sample);
    return maxHalfplaneLin(net_points.begin(), net_points.end(), sample_points.begin(), sample_points.end());
}

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

    py::class_<pyscan::Grid<int>>("BiGrid", py::no_init)
            .def("totalRedWeight", &pyscan::Grid<int>::totalRedWeight)
            .def("totalBlueWeight", &pyscan::Grid<int>::totalBlueWeight)
            .def("redWeight", &pyscan::Grid<int>::redWeight)
            .def("blueWeight", &pyscan::Grid<int>::blueWeight)
            .def("xCoord", &pyscan::Grid<int>::xCoord)
            .def("yCoord", &pyscan::Grid<int>::yCoord)
            .def("toRectangle", &pyscan::Grid<int>::toRectangle)
            .def("size", &pyscan::Grid<int>::size);

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

    py::class_<pyscan::Point<int>>("Point", py::init<double, double, int, int>())
            .def("getWeight", &pyscan::Point<int>::getWeight)
            .def("getRedWeight", &pyscan::Point<int>::getRedWeight)
            .def("getBlueWeight", &pyscan::Point<int>::getBlueWeight)
            .def("getX", &pyscan::Point<int>::getX)
            .def("getY", &pyscan::Point<int>::getY)
            .def("__str__", &pyscan::Point<int>::toString)
            .def("__repr__", &pyscan::Point<int>::toString);

    py::class_<BloomFilter>("BloomFilter", py::init<int, double>())
            .def("insert", &BloomFilter::insert)
            .def("mightBePresent", &BloomFilter::mightBePresent);

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

    py::class_<pyscan::Halfplane>("Halfplane", py::init<double, double, double>())
            .def("getSlope", &pyscan::Halfplane::getSlope)
            .def("getIntersect", &pyscan::Halfplane::getIntersect)
            .def("fValue", &pyscan::Halfplane::fValue);

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
    py::def("maxSubgridKullSlow", &pyscan::maxSubgridKullSlow<int>);
    py::def("maxSubgridLinearSlow", &pyscan::maxSubgridLinearSlow<int>);
    py::def("maxSubgridLinearSimple", &pyscan::maxSubgridLinearSimple<int>);
    py::def("maxSubgridLinear", &pyscan::maxSubgridLinearG<int>);

    py::def("maxHalfPlaneStat", &pyscan::maxHalfplaneStat);
    py::def("maxHalfPlaneLin", &pyscan::maxHalfplaneLin);
    py::def("maxHalfPlaneGamma", &pyscan::maxHalfplaneGamma);
}

