//
// Created by mmath on 7/7/17.
//

#include <functional>
#include <iostream>


#include <boost/iterator_adaptors.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>



#include "RectangleScan.hpp"
#include "HalfplaneScan.hpp"
#include "DiskScan.hpp"
#include "DiskScan2.hpp"
#include "FunctionApprox.hpp"
#include "TrajectoryScan.hpp"

namespace py = boost::python;

template <class K, class V>
boost::python::dict toPythonDict(std::unordered_map<K, V> const& map) {
    boost::python::dict dictionary;
    for (auto iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

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

template<class T>
py::list std_vector_to_py_list(const std::vector<T>& v) {
    py::object get_iter = py::iterator<std::vector<T> >();
    py::object iter = get_iter(v);
    py::list l(iter);
    return l;
}


namespace pyscan {

    auto sized_region(double size) -> std::function<double(double, double)> {
        /*
         * Useful for finding a region of a certain size.
         */

        return [size] (double m, double b) {
            //std::cout << 1 - fabs(m - size) << " " << size << std::endl;
            return 1 - fabs(m - size);
        };
    }

    auto toLPt(const py::object& pts, 
            const py::object& weights,
            const py::object& labels) -> pyscan::lpoint_list {

        auto bp = py::stl_input_iterator<Pt2>(pts);
        auto ep = py::stl_input_iterator<Pt2>();
        auto bw = py::stl_input_iterator<double>(weights);
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

    Subgrid maxSubgridLin(Grid const& grid, double eps, std::function<double(double, double)> const& f) {
      return maxSubgridLinearSimple(grid, eps, f);
    }

    Subgrid maxSubgridLinTheory(Grid const& grid, double eps, std::function<double(double, double)> const& f){
        return maxSubgridLinearTheory(grid, eps, f);
    }

    Subgrid maxSubgridSlow(Grid const &grid, std::function<double(double, double)> const& f) {
        return maxSubgridNonLinear(grid, f);
    }

    py::tuple maxHalfspace(const py::object& net,
                        const py::object& sampleM,
                        const py::object& weightM,
                        const py::object& sampleB,
                        const py::object& weightB,
                        std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto red = to_std_vector<Pt2>(sampleM);
        auto red_w = to_std_vector<double>(weightM);
        auto blue = to_std_vector<Pt2>(sampleB);
        auto blue_w = to_std_vector<double>(weightB);
        Pt2 d1;
        double d1value;
        std::tie(d1, d1value) = max_halfplane(net_points,
                red, red_w, blue, blue_w, f);
        return py::make_tuple(d1, d1value);
    }

    py::tuple maxHalfspaceLabels(const py::object& net,
                           const py::object& sampleM,
                           const py::object& weightM,
                           const py::object& labelM,
                           const py::object& sampleB,
                           const py::object& weightB,
                           const py::object& labelB,
                           std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto red = to_std_vector<Pt2>(sampleM);
        auto red_w = to_std_vector<double>(weightM);
        auto red_l = to_std_vector<long>(labelM);
        auto blue = to_std_vector<Pt2>(sampleB);
        auto blue_w = to_std_vector<double>(weightB);
        auto blue_l = to_std_vector<long>(labelB);
        Pt2 d1;
        double d1value;
        std::tie(d1, d1value) = max_halfplane_labeled(net_points,
                                              red, red_w, red_l,
                                                      blue, blue_w, blue_l, f);
        return py::make_tuple(d1, d1value);
    }
    py::tuple maxDisk(const py::object& net,
            const py::object& sampleM, 
            const py::object& weightM,
            const py::object& sampleB,
            const py::object& weightB,
            std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto sample_p_M = toWPt(sampleM, weightM);
        auto sample_p_B = toWPt(sampleB, weightB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = diskScan(net_points, sample_p_M, 
            sample_p_B, 
            f);
        return py::make_tuple(d1, d1value);
    }

//    py::tuple maxDiskScale(const py::object& net,
//                           const py::object& sampleM,
//                           const py::object& sampleB,
//                           int r,
//                           std::function<double(double, double)> const& f) {
//        auto net_points = to_std_vector<pyscan::Point<>>(net);
//        auto m_sample = to_std_vector<pyscan::WPoint<>>(sampleM);
//        auto b_sample = to_std_vector<pyscan::WPoint<>>(sampleB);
//        Disk d1;
//        double d1value;
//        std::tie(d1, d1value) = disk_scan_scale(net_points,
//                                                m_sample,
//                                                b_sample,
//                                                r,
//                                                f);
//        return py::make_tuple(d1, d1value);
//    }

    py::tuple maxDiskCached(const py::object& net,
                           const py::object& sampleM,
                           const py::object& sampleB,
                           std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto m_sample = to_std_vector<pyscan::WPoint<>>(sampleM);
        auto b_sample = to_std_vector<pyscan::WPoint<>>(sampleB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = cached_disk_scan(net_points,
                                                m_sample,
                                                b_sample,
                                                f);
        return py::make_tuple(d1, d1value);
    }
    py::tuple maxDiskLCached(const py::object& net,
                            const py::object& sampleM,
                            const py::object& sampleB,
                            std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto m_sample = to_std_vector<pyscan::LPoint<>>(sampleM);
        auto b_sample = to_std_vector<pyscan::LPoint<>>(sampleB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = cached_disk_scan(net_points,
                                                 m_sample,
                                                 b_sample,
                                                 f);
        return py::make_tuple(d1, d1value);
    }
//    py::tuple maxDiskScaleLabel(const py::object& net,
//                           const py::object& sampleM,
//                           const py::object& sampleB,
//                           int r,
//                           std::function<double(double, double)> const& f) {
//        auto net_points = to_std_vector<pyscan::Point<>>(net);
//        auto m_sample = to_std_vector<pyscan::LPoint<>>(sampleM);
//        auto b_sample = to_std_vector<pyscan::LPoint<>>(sampleB);
//        Disk d1;
//        double d1value;
//        std::tie(d1, d1value) = label_disk_scan_scale(net_points,
//                                                m_sample,
//                                                b_sample,
//                                                r,
//                                                f);
//        return py::make_tuple(d1, d1value);
//    }

    py::tuple maxDiskLabels(const py::object& net,
            const py::object& sampleM, 
            const py::object& weightM,
            const py::object& labelM,
            const py::object& sampleB,
            const py::object& weightB,
            const py::object& labelB,
            std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto lpointsM = toLPt(sampleM, weightM, labelM);
        auto lpointsB = toLPt(sampleB, weightB, labelB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = diskScanLabels(net_points, lpointsM, lpointsB, f);
        return py::make_tuple(d1, d1value);
    }

      py::tuple maxHalfplane(const py::object& net, 
            const py::object& sampleM, 
            const py::object& weightM,
            const py::object& sampleB,
            const py::object& weightB,
            std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<pyscan::Point<>>(net);
        auto sample_p_M = to_std_vector<Pt2>(sampleM);
        auto sample_p_B = to_std_vector<Pt2>(sampleB);
        auto weight_p_M = to_std_vector<double>(weightM);
        auto weight_p_B = to_std_vector<double>(weightB);
        Point<> d1;
        double d1value;
        std::tie(d1, d1value) = max_halfplane(net_points, 
            sample_p_M,
            weight_p_M,
            sample_p_B,
            weight_p_B, 
            f);
        return py::make_tuple(d1, d1value);
    }

    py::tuple maxHalfplaneLabeled(const py::object& net, 
            const py::object& sampleM, 
            const py::object& weightM,
            const py::object& labelM,
            const py::object& sampleB,
            const py::object& weightB,
            const py::object& labelB,
            std::function<double(double, double)> const& f) {
        auto net_points = to_std_vector<Pt2>(net);
        auto lpointsM = toLPt(sampleM, weightM, labelM);
        auto lpointsB = toLPt(sampleB, weightB, labelB);
        Disk d1;
        double d1value;
        std::tie(d1, d1value) = diskScanLabels(net_points, lpointsM, lpointsB, f);
        return py::make_tuple(d1, d1value);
    }

    py::object approx_hull(double eps, std::function<double(pyscan::Vec2)> const& phi, const py::object& traj_pts) {
        auto pts = to_std_vector<Pt2>(traj_pts);
        auto max_f = [&] (Vec2 direction) {
            double max_dir = -std::numeric_limits<double>::infinity();
            Pt2 curr_pt {0.0, 0.0, 0.0};
            for (auto& pt : pts) {
                double curr_dir = direction[0] * pyscan::getX(pt) + direction[1] * pyscan::getY(pt);
                if (max_dir < curr_dir) {
                    max_dir = curr_dir;
                    curr_pt = pt;
                }
            }
            return Vec2{pyscan::getX(curr_pt), pyscan::getY(curr_pt)};
        };
        std::vector<pyscan::Point<>> core_set_pts;
        {
            auto vecs = eps_core_set(eps, max_f);
            for (auto &v :vecs) {
                core_set_pts.push_back(pyscan::Point<>(v[0], v[1], 1.0));
            }
        }
        return std_vector_to_py_list(core_set_pts);
    }


    py::object approx_hull3(double eps, std::function<double(pyscan::Vec3)> const& phi, const py::object& traj_pts) {
        /*
         * Finish this.
         */
        auto pts = to_std_vector<Pt3>(traj_pts);
        auto max_f = [&] (Vec3 direction) {
            double max_dir = -std::numeric_limits<double>::infinity();
            Pt3 curr_pt {0.0, 0.0, 0.0, 0.0};
            for (auto& pt : pts) {
                double curr_dir = direction[0] * pyscan::getX(pt)
                        + direction[1] * pyscan::getY(pt)
                        + direction[2] * pyscan::getZ(pt);
                if (max_dir < curr_dir) {
                    max_dir = curr_dir;
                    curr_pt = pt;
                }
            }
            return Vec3{pyscan::getX(curr_pt), pyscan::getY(curr_pt), pyscan::getZ(curr_pt)};
        };
        std::vector<pyscan::Point<3>> core_set_pts;
        {
            auto vecs = eps_core_set3(eps, max_f);
            for (auto &v :vecs) {
                core_set_pts.push_back(pyscan::Point<3>(v[0], v[1], v[2], 1.0));
            }
        }
        return std_vector_to_py_list(core_set_pts);
    }

    double evaluate(std::function<double(double, double)> const& f, double m, double b) {
        return f(m, b);
    }



    pyscan::traj_set create_traj_set(const py::object& net) {
        /*
         * Expects a list of lists.
         */
        auto beg_net = py::stl_input_iterator<py::list>(net);
        auto end_net = py::stl_input_iterator<py::list>();
        pyscan::point_list net_traj_pts;
        std::vector<size_t> net_offsets;
        size_t next_offset = 0;
        for (auto b = beg_net; b != end_net; b++) {
            auto curr_traj = py::stl_input_iterator<Point<>>(*b);
            auto end_traj = py::stl_input_iterator<Point<>>();
            for (auto b_traj = curr_traj; b_traj != end_traj; b_traj++) {
                next_offset += 1;
                net_traj_pts.push_back(*b_traj);
            }
            net_offsets.push_back(next_offset);
        }
        return {net_traj_pts, net_offsets};
    }

    pyscan::wtraj_set create_wtraj_set(const py::object& pts, const py::object& weights) {
        pyscan::traj_set trajectories = create_traj_set(pts);
        auto cv_weights = to_std_vector<double>(weights);
        return {trajectories.traj_pts, trajectories.offsets, cv_weights};
    }

    py::tuple traj_disk_scan_py(const py::object& net,
                            const py::object& sampleM,
                            const py::object& weightM,
                            const py::object& sampleB,
                            const py::object& weightB,
                            double alpha,
                            double min_r,
                            std::function<double(double, double)> const& f) {

        auto net_set = create_traj_set(net);
        auto m_set = create_wtraj_set(sampleM, weightM);
        auto b_set = create_wtraj_set(sampleB, weightB);

        Disk d1;
        double d1value;
        std::tie(d1, d1value) = pyscan::traj_disk_scan(net_set, m_set, b_set, alpha, min_r, f);
        return py::make_tuple(d1, d1value);
    }

    py::dict grid_traj_py(py::list const& traject, double chord_l) {
        auto points = to_std_vector<Point<>>(traject);
        auto mapped_pts = pyscan::grid_traj(points.begin(), points.end(), chord_l);
        return toPythonDict(mapped_pts);
    }

    py::dict approx_traj_cells_py(py::list const& traject, double chord_l, double eps) {
        auto points = to_std_vector<Point<>>(traject);
        auto mapped_pts = pyscan::approximate_traj_cells(points.begin(), points.end(),
                chord_l, eps);
        return toPythonDict(mapped_pts);
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
/// @brief Type that allows for registration of conversions from
///        python iterable types.
struct iterable_converter
{
    /// @note Registers converter from a python interable type to the
    ///       provided type.
    template <typename Container>
    iterable_converter&
    from_python()
    {
        boost::python::converter::registry::push_back(
                &iterable_converter::convertible,
                &iterable_converter::construct<Container>,
                boost::python::type_id<Container>());

        // Support chaining.
        return *this;
    }

    /// @brief Check if PyObject is iterable.
    static void* convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    /// @brief Convert iterable PyObject to C++ container type.
    ///
    /// Container Concept requirements:
    ///
    ///   * Container::value_type is CopyConstructable.
    ///   * Container can be constructed and populated with two iterators.
    ///     I.e. Container(begin, end)
    template <typename Container>
    static void construct(
            PyObject* object,
            boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace python = boost::python;
        // Object is a borrowed reference, so create a handle indicting it is
        // borrowed for proper reference counting.
        python::handle<> handle(python::borrowed(object));

        // Obtain a handle to the memory block that the converter has allocated
        // for the C++ type.
        typedef python::converter::rvalue_from_python_storage<Container>
                storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        typedef python::stl_input_iterator<typename Container::value_type>
                iterator;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        new (storage) Container(
                iterator(python::object(handle)), // begin
                iterator());                      // end
        data->convertible = storage;
    }
};

struct pypoint_converter {


    /// @note Registers converter from a python iterable type to the
    ///       provided type.
    pypoint_converter& from_python()
    {
        boost::python::converter::registry::push_back(
                &pypoint_converter::convertible,
                &pypoint_converter::construct,
                boost::python::type_id<pyscan::Point<>>());
        return *this;
    }

    /// @brief Check if PyObject is a string.
    static void* convertible(PyObject* object) {
        if (PyTuple_Check(object) && PyTuple_Size(object) == 2)  {
            return PyFloat_Check(PyTuple_GetItem(object, 0)) &&
            PyFloat_Check(PyTuple_GetItem(object, 1)) ? object : NULL;
        }
        return NULL;
    }

    /// @brief Convert PyString to Container.
    ///
    /// Container Concept requirements:
    ///
    ///   * Container::value_type is CopyConstructable from char.
    ///   * Container can be constructed and populated with two iterators.
    ///     I.e. Container(begin, end)
    static void construct(
            PyObject* object,
            boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace python = boost::python;
        // Object is a borrowed reference, so create a handle indicting it is
        // borrowed for proper reference counting.
        python::handle<> handle(python::borrowed(object));

        // Obtain a handle to the memory block that the converter has allocated
        // for the C++ type.
        typedef python::converter::rvalue_from_python_storage<pyscan::Point<>>
                storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        double x = PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 0));
        double y = PyFloat_AsDouble(PyTuple_GetItem(object, 1));
        data->convertible = new (storage) pyscan::Point<>(x, y, 1.0);
    }
};


template<typename T1, typename T2>
struct tuple_to_python_tuple {

    static PyObject* convert(const std::tuple<T1, T2> &s) {
        return boost::python::incref(py::make_tuple(std::get<0>(s), std::get<1>(s)).ptr());
    }
};

void test1(std::vector<pyscan::Point<>> values)
{
    for (auto&& value: values)
        std::cout << value << std::endl;
}


BOOST_PYTHON_MODULE(libpyscan) {
    using namespace py;
    /*
    py::class_<pyscan::MaximumIntervals>("MaximumIntervals", py::init<std::size_t>())
            .def("mergeZeros")
    */

    //Should convert tuples directly pyscan points.
    pypoint_converter().from_python();

    // Register interable conversions.
    iterable_converter()
            // Build-in type.
            .from_python<std::vector<double> >()
                    // Each dimension needs to be convertable.
            .from_python<std::vector<pyscan::Point<> >>()
            .from_python<std::vector<std::vector<pyscan::Point<> > > >();

    to_python_converter<std::tuple<pyscan::Disk, double>, tuple_to_python_tuple<pyscan::Disk, double>>();
    to_python_converter<std::tuple<pyscan::Pt2, double>, tuple_to_python_tuple<pyscan::Pt2, double>>();

    py::def("test1", &test1);
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

    py::class_<pyscan::Pt2>("Point", py::init<double, double, double>())
            .def("approx_eq", &pyscan::Point<>::approx_eq)
            .def("__getitem__", &pyscan::Point<>::operator[])
            .def("above", &pyscan::Point<>::above)
            .def("above_closed", &pyscan::Point<>::above_closed)
            .def("__str__", &pyscan::Point<>::str)
            .def("__repr__", &pyscan::Point<>::str)
            .def("__eq__", &pyscan::Point<>::operator==)
            .def("getX", &pyscan::getX<2>)
            .def("getY", &pyscan::getY<2>);

    py::class_<pyscan::WPoint<2>, py::bases<pyscan::Pt2>>("WPoint", py::init<double, double, double, double>())
            .def("get_weight", &pyscan::WPoint<2>::get_weight);

    py::class_<pyscan::LPoint<2>, py::bases<pyscan::WPoint<2>>>("LPoint", py::init<size_t, double, double, double, double>())
            .def("get_label", &pyscan::LPoint<2>::get_label);

    py::class_<std::function<double(double, double)> >("CFunction", py::no_init);

    py::scope().attr("KULLDORF") = std::function<double(double, double)>(
        [&](double m, double b) {
            return pyscan::kulldorff(m, b, .0001);
    });

    py::scope().attr("DISC") = std::function<double(double, double)>(
        [&](double m, double b) {
            return fabs(m - b);
    });

    py::scope().attr("RKULLDORF") = std::function<double(double, double)>(
        [&](double m, double b) {
            return pyscan::regularized_kulldorff(m, b, .0001);
    });

    py::def("evaluate", &pyscan::evaluate);
    py::def("size_region", &pyscan::sized_region);

    //py::def("dot", &pyscan::dot<2ul>);
    py::def("intersection", &pyscan::intersection);
    py::def("correct_orientation", &pyscan::correct_orientation);


    py::class_<pyscan::Disk>("Disk", py::init<double, double, double>())
            .add_property("a", &pyscan::Disk::getA, &pyscan::Disk::setA)
            .add_property("b", &pyscan::Disk::getB, &pyscan::Disk::setB)
            .add_property("radius", &pyscan::Disk::getR, &pyscan::Disk::setR)
            .def("contains", &pyscan::Disk::contains);


    py::def("make_grid", makeGrid);
    //py::def("makeNetGrid", makeNetGrid);
    //py::def("makeSampleGrid", makeSampleGrid);

    py::def("max_subgrid_slow", &pyscan::maxSubgridSlow);
    py::def("max_subgrid_linear_simple", &pyscan::maxSubgridLin);
    py::def("max_subgrid_linear", &pyscan::maxSubgridLinTheory);

    //py::def("maxDiskLabels", &maxDiskLabelsI);
    py::def("max_halfplane", &pyscan::max_halfplane);
    py::def("max_halfplane_labels", &pyscan::maxHalfspaceLabels);
    py::def("max_disk", &pyscan::maxDisk);
    py::def("max_disk_labels", &pyscan::maxDiskLabels);

    py::def("max_traj_disk", &pyscan::traj_disk_scan_py);
    py::def("grid_traj", &pyscan::grid_traj_py);
    py::def("approx_traj_cells", &pyscan::approx_traj_cells_py);

    //   py::def("max_disk_scale_labels", &pyscan::maxDiskScaleLabel);

    py::def("max_disk_cached", &pyscan::maxDiskCached);
    py::def("max_disk_label_cached", &pyscan::maxDiskLCached);

    //py::def("maxRectLabels", &maxRectLabelsI);
    //py::def("maxRectLabels", &maxRectLabelsD);
    py::def("approximate_hull", pyscan::approx_hull);

    py::def("approximate_hull3", pyscan::approx_hull3);

}
