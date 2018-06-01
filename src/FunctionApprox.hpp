//
// Created by mmath on 7/7/17.
//

#ifndef PYSCAN_FUNCTIONAPPROX_HPP
#define PYSCAN_FUNCTIONAPPROX_HPP
#include <functional>
#include <cmath>
#include <array>

namespace pyscan {

	using Vec2 = std::array<double, 2>;

	inline double dot(Vec2 const& v1, Vec2 const& v2) {
      return v1[0] * v2[0] + v1[1] * v2[1];
    }

    inline double mag(Vec2 const& v1) {
    	return sqrt(dot(v1, v1));
    }

    inline Vec2 operator+(Vec2 const& v1, Vec2 const& v2) {
      return {v1[0] + v1[0], v2[1] + v2[1]};
    }

    inline Vec2 operator*(Vec2 const& v2, double val) {
    	return {val * v2[0], val * v2[1]};
    }

    inline Vec2 operator/(Vec2 const& v2, double val) {
    	return {v2[0] / val, v2[1] / val};
    }

	double approximateHull(double eps,
							Vec2 const& cc, Vec2 const& cl,
							std::function<double(Vec2)> phi, //function to maximize
							std::function<Vec2(Vec2)> lineMaxF);

	double approximateHull(double eps,
	        	std::function<double(Vec2)> phi, //function to maximize
	            std::function<Vec2(Vec2)> lineMaxF);
}
#endif //PYSCAN_FUNCTIONAPPROX_HPP
