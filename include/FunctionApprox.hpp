//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_FUNCTIONAPPROX_HPP
#define PYSCAN_FUNCTIONAPPROX_HPP
#include <functional>
#include <cmath>
#include <vector>
#include <array>

#include "Point.hpp"

namespace pyscan {

    template<size_t dim>
	using Vec = std::array<double, dim>;

    using Vec2 = Vec<2>;
    using Vec3 = Vec<3>;

    template <size_t dim>
	inline double dot(Vec<dim> const& v1, Vec<dim> const& v2) {
        double accum = 0;
        for (size_t i = 0; i < dim; i++) {
            accum += v1[i] * v2[i];
        }
      return accum;
    }

    template <size_t dim>
    inline double mag(Vec<dim> const& v1) {
    	return sqrt(dot(v1, v1));
    }

    template <size_t dim>
    inline Vec<dim> operator+(Vec<dim> const& v1, Vec<dim> const& v2) {
        Vec<dim> v_out;
        for (size_t i = 0; i < dim; i++) {
            v_out[i] = v1[i] + v2[i];
        }
        return v_out;
    }

    template <size_t dim>
    inline Vec<dim> operator*(Vec<dim> const& v2, double val) {
        Vec<dim> v_out;
        for (size_t i = 0; i < dim; i++) {
            v_out[i] = v2[i] * val;
        }
        return v_out;
    }

    template<size_t dim>
    inline Vec<dim> operator/(Vec<dim> const& v2, double val) {
        Vec<dim> v_out;
        for (size_t i = 0; i < dim; i++) {
            v_out[i] = v2[i] / val;
        }
        return v_out;
    }


	double approximateHull(double eps,
							Vec2 const& cc, Vec2 const& cl,
							std::function<double(Vec2)> phi, //function to maximize
							std::function<Vec2(Vec2)> lineMaxF);

	double approximateHull(double eps,
	        	std::function<double(Vec2)> phi, //function to maximize
	            std::function<Vec2(Vec2)> lineMaxF);

	std::vector<Vec2> eps_core_set(double eps,
								   std::function<Vec2(Vec2)> lineMaxF);

    point_list_t approx_hull(point_list_t const &pts, double eps);
}
#endif //PYSCAN_FUNCTIONAPPROX_HPP
