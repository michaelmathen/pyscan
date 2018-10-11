//
// Created by mmath on 7/7/17.
//

#include <functional>
#include <iostream>


#include <boost/iterator_adaptors.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>



#include "RectangleScan.hpp"
#include "HalfSpaceScan.hpp"
#include "DiskScan.hpp"
#include "DiskScan2.hpp"
#include "FunctionApprox.hpp"
#include "TrajectoryScan.hpp"
#include "TrajectoryCoreSet.hpp"

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

auto toPointList(py::object const& el) -> pyscan::point_list_t {
    pyscan::point_list_t l;
    for (auto beg = py::stl_input_iterator<pyscan::Point<>>(el);
         beg != py::stl_input_iterator<pyscan::Point<>>(); beg++) {
        l.push_back(*beg);
    }
    return l;
}

template<class T>
py::list std_vector_to_py_list(const std::vector<T>& v) {
    py::list l;
    for (auto it = v.begin(); it != v.end(); it++) {
        l.append(*it);
    }
    return l;
}


namespace pyscan {

    auto sized_region(double size) -> std::function<double(double, double)> {
        /*
         * Useful for finding a region of a certain size.
         */
        return [size] (double m, double b) {
            (void)b;
            return 1 - fabs(m - size);
        };
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


    double evaluate(std::function<double(double, double)> const& f, double m, double b) {
        return f(m, b);
    }

};



pyscan::Grid makeGrid(const py::object& sample_r, const py::object& weight_r, const py::object& sample_b,
                      const py::object& weight_b, size_t r) {
    auto points_r = to_std_vector<pyscan::Point<>>(sample_r);
    auto points_b = to_std_vector<pyscan::Point<>>(sample_b);
    auto v_weight_r = to_std_vector<double>(weight_r);
    auto v_weight_b = to_std_vector<double>(weight_b);
    return pyscan::Grid(r, points_r, v_weight_r, points_b, v_weight_b);
}

struct iterable_converter {
    template <typename Container>
    iterable_converter&
    from_python()
    {
        boost::python::converter::registry::push_back(
                &iterable_converter::convertible,
                &iterable_converter::construct<Container>,
                boost::python::type_id<Container>());
        return *this;
    }

    static void* convertible(PyObject* object)
    {
        return PyObject_GetIter(object) ? object : NULL;
    }

    template <typename Container>
    static void construct(
            PyObject* object,
            boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        namespace python = boost::python;
        python::handle<> handle(python::borrowed(object));

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


template<int dim>
struct pypoint_converter {

    static_assert(dim == 3 || dim == 2, "Converter not implemented for dimensions other than 2 or 3.");

    pypoint_converter& from_python() {
        boost::python::converter::registry::push_back(
                &pypoint_converter::convertible,
                &pypoint_converter::construct,
                boost::python::type_id<pyscan::Point<dim>>());
        return *this;
    }

    /// @brief Check if PyObject is a double tuple.
    static void* convertible(PyObject* object) {
        if (PyTuple_Check(object) && PyTuple_Size(object) == dim)  {
            for (int i = 0; i < dim; i++) {
               if (!PyFloat_Check(PyTuple_GetItem(object, 0))) {
                   return NULL;
               }
            }
            return object;
        }
        return NULL;
    }

    static void construct( PyObject* object, boost::python::converter::rvalue_from_python_stage1_data* data) {
        namespace python = boost::python;
        python::handle<> handle(python::borrowed(object));
        typedef python::converter::rvalue_from_python_storage<pyscan::Point<dim>>
                storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        if (dim == 2) {
            data->convertible = new (storage) pyscan::Point<>(PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 0)),
                                                                PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 1)), 1.0);
        } else {
            data->convertible = new (storage) pyscan::Point<3>(PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 0)),
                                                              PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 1)),
                                                               PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 2)), 1.0);
        }
    }
};


struct pytrajectory_converter {

    pytrajectory_converter& from_python() {
        boost::python::converter::registry::push_back(
                &pytrajectory_converter::convertible,
                &pytrajectory_converter::construct,
                boost::python::type_id<pyscan::wtrajectory_t>());
        return *this;
    }

    /// @brief Check if PyObject is a double tuple.
    static void* convertible(PyObject* object) {
        if (PyTuple_Check(object) && PyTuple_Size(object) == 2) {
            if (PyFloat_Check(PyTuple_GetItem(object, 0)) && PyIter_Check(PyTuple_GetItem(object, 1))) {
                return object;
            }
        }
        return NULL;
    }

    static void construct( PyObject* object, boost::python::converter::rvalue_from_python_stage1_data* data) {
        namespace python = boost::python;
        python::handle<> handle(python::borrowed(object));
        typedef python::converter::rvalue_from_python_storage<pyscan::wtrajectory_t>
                storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.  The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.

        data->convertible = new (storage) pyscan::WTrajectory(PyFloat_AS_DOUBLE(PyTuple_GetItem(object, 0)),
                py::extract<pyscan::point_list_t>(PyTuple_GetItem(object, 1)));
    }
};



/*
 * This converts two argument c++ tuples into python tuples automatically when they are returned from the c++ code.
 * s -- a two element tuple of type t1 and t2
 * returns -- a two element python tuple of type t1 and t2.
 */
template<typename T1, typename T2>
struct tuple_to_python_tuple {
    static PyObject* convert(const std::tuple<T1, T2> &s) {
        return boost::python::incref(py::make_tuple(std::get<0>(s), std::get<1>(s)).ptr());
    }
};


/*
 * This converts returned vectors automatically into python lists when the are returned from the c++ code.
 * s -- a vector of type T
 * returns -- A python list containing the elements of s.
 */
template <typename T>
struct vector_to_python_list {
    static PyObject* convert(const std::vector<T> &s) {
        auto new_list = std_vector_to_py_list(s);
        return boost::python::incref(new_list.ptr());
    }
};

BOOST_PYTHON_MODULE(libpyscan) {
    using namespace py;
    /*
    py::class_<pyscan::MaximumIntervals>("MaximumIntervals", py::init<std::size_t>())
            .def("mergeZeros")
    */

    //Should convert tuples directly pyscan points.
    pypoint_converter<2>().from_python();
    pypoint_converter<3>().from_python();

    pytrajectory_converter().from_python();

    // Register interable conversions.
    iterable_converter()
            // Build-in type.
            .from_python<std::vector<double> >()
                    // Each dimension needs to be convertable.
            .from_python<std::vector<pyscan::Point<> >>()
            .from_python<std::vector<pyscan::WPoint<>>>()
            .from_python<std::vector<pyscan::LPoint<>>>()
            .from_python<std::vector<std::vector<pyscan::Point<> > > >()
            .from_python<std::vector<std::vector<pyscan::WPoint<> > > >();

    to_python_converter<std::tuple<pyscan::Disk, double>, tuple_to_python_tuple<pyscan::Disk, double>>();
    to_python_converter<std::tuple<pyscan::pt2_t, double>, tuple_to_python_tuple<pyscan::pt2_t, double>>();

    to_python_converter<std::vector<double>, vector_to_python_list<double>>();
    to_python_converter<std::vector<pyscan::Point<2>>, vector_to_python_list<pyscan::Point<2>>>();
    to_python_converter<std::vector<pyscan::Point<3>>, vector_to_python_list<pyscan::Point<3>>>();

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

    py::class_<pyscan::pt2_t>("Point", py::init<double, double, double>())
            .def("approx_eq", &pyscan::Point<>::approx_eq)
            .def("__getitem__", &pyscan::Point<>::operator())
            .def("above", &pyscan::Point<>::above)
            .def("above_closed", &pyscan::Point<>::above_closed)
            .def("__str__", &pyscan::Point<>::str)
            .def("__repr__", &pyscan::Point<>::str)
            .def("__eq__", &pyscan::Point<>::operator==);

    py::class_<pyscan::WPoint<2>, py::bases<pyscan::pt2_t>>("WPoint", py::init<double, double, double, double>())
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
            .add_property("a", &pyscan::Disk::getA)
            .add_property("b", &pyscan::Disk::getB)
            .add_property("radius", &pyscan::Disk::getR)
            .def("contains", &pyscan::Disk::contains);


    py::def("make_grid", makeGrid);
    //py::def("makeNetGrid", makeNetGrid);
    //py::def("makeSampleGrid", makeSampleGrid);

    py::def("max_subgrid_slow", &pyscan::maxSubgridSlow);
    py::def("max_subgrid_linear_simple", &pyscan::maxSubgridLin);
    py::def("max_subgrid_linear", &pyscan::maxSubgridLinTheory);


    //Max Halfspace codes
    py::def("max_halfplane", &pyscan::MaxHalfPlane);
    py::def("max_halfplane_labels", &pyscan::MaxHalfPlaneLabeled);
    py::def("max_halfspace", &pyscan::MaxHalfSpace);
    py::def("max_halfspace_labels", &pyscan::MaxHalfSpaceLabeled);


    py::def("max_disk", &pyscan::disk_scan);
    py::def("max_disk_labels", &pyscan::disk_scan_labels);


    //   py::def("max_disk_scale_labels", &pyscan::maxDiskScaleLabel);

    py::def("max_disk_cached", &pyscan::cached_disk_scan<pyscan::wpoint_list_t >);
    py::def("max_disk_label_cached", &pyscan::cached_disk_scan<pyscan::lpoint_list_t >);

    py::def("max_rect_labels", &pyscan::max_rect_labels);

    ////////////////////////////////////////////////////////////////////
    //TrajectoryScan.hpp wrappers///////////////////////////////////////
    ////////////////////////////////////////////////////////////////////

    py::class_<pyscan::wtrajectory_t>("WTrajectory", py::init<double, pyscan::point_list_t>())
            .def("get_weight", &pyscan::wtrajectory_t::get_weight);

    py::def("max_disk_traj_grid", &pyscan::max_disk_traj_grid);


    //This grids the trajectory and assigns a single point to each cell.
    py::def("grid_kernel", &pyscan::approx_traj_grid);

    //This grids the trajectory and creates an alpha hull in each one.
    py::def("grid_direc_kernel", &pyscan::approx_traj_kernel_grid);

    //This is for 2d eps-kernel useful for halfspaces.
    py::def("approximate_hull", pyscan::approx_hull);

    //This is a 3d eps-kernel for disks.
    py::def("lifting_kernel", &pyscan::lifting_coreset);
}
