//
// Created by mmath on 2/12/18.
//
#include <tuple>
#include <algorithm>
#include <cstddef>
#include "Cuttings.hpp"

namespace pyscan {

    std::tuple<double, double, double>  cross_product(double a1, double b1, double c1,
                                                      double a2, double b2, double c2) {
        return std::make_tuple(b1 * c2 - b2 * c1, a2 * c1 - a1 * c2, a1 * b2 - a2 * b1);
    }

	ProjObject ProjObject::intersection(ProjObject const& other_line) {
        auto args = cross_product(this->a, this->b, this->c, other_line.a, other_line.b, other_line.c);
		return ProjObject(std::get<0>(args), std::get<1>(args), std::get<2>(args));
	}

    bool ProjObject::below(ProjObject const& obj1){
        return obj1.a * this->a + obj1.b * this->b + obj1.c * this->c < 0;
    }

    bool ProjObject::below_closed(ProjObject const& obj1){
        return obj1.a * this->a + obj1.b * this->b + obj1.c * this->c <= 0;
    }

    bool ProjObject::order(ProjObject const& l1, ProjObject const& l2){
        /*
         * Orders the intersection of two lines with the current line via there x coordinates.
         * Or computes two lines going through the same point passing through l1 and l2 and orders
         * them by their slope.
         */
        auto p1 = this->intersection(l1);
        auto p2 = this->intersection(l2);
        return p1->x_order(p2);
    }

    bool ProjObject::x_order(ProjObject const& p2) {
        /*
         * Uses the x_coordinate to figure out which point is lower.
         * (or tests to see if a line has a steeper slope).
         */
        if (0 <= this->c * p2.c ) {
            return this->a * p2.c < p2.a * this->c;
        } else {
            return this->a * p2.c > p2.a * this->c;
        }
    }

    void Cell::insert_points(const std::vector<std::tuple<ProjObject, int>> &pts) {

    }

    void Cell::insert_lines(const std::vector<ProjObject> &lines) {

    }

    bool Triangle::intersects(const ProjObject &l1) {
        return false;
    }

    std::vector<Cell> triangle_cutting(const std::vector<ProjObject> &lines, int r) {

    	using inters_t = std::tuple<ProjObject*, ProjObject*>;
        std::vector<ProjObject> internal_lines = lines;

        //Add two lines at infinity. The return points at opposite infinities.
        internal_lines.push_back(ProjObject(0, 0, 1));
        internal_lines.push_back(ProjObject(0, 0, -1));
        std::vector<inters_t> x_intercepts;

        for (int i = 0; i < internal_lines.size() - 1; i++ ) {
			for (int j = i + 1; j < internal_lines.size(); j++) {
				x_intercepts.emplace_back(&internal_lines[i], &internal_lines[j]);
			}
        }


        std::sort(x_intercepts.begin(), x_intercepts.end(), [&](inters_t const& isp1, inters_t const& isp2) {
            auto inters_p1 = std::get<0>(isp1)->intersection(*std::get<1>(isp1));
            auto inters_p2 = std::get<0>(isp2)->intersection(*std::get<1>(isp2));
        	return inters_p1.x_order(inters_p2);
        });

        std::vector<Triangle> curr_triangles;
        std::vector<ProjObject*> first_lines(internal_lines.size(), nullptr);
        std::vector<ProjObject*> second_lines(internal_lines.size(), nullptr);
        for (auto& intersect : x_intercepts) {

        }

    }
}