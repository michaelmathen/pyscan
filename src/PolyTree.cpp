//
// Created by mmath on 2/12/18.
//
#include <tuple>
#include <algorithm>
#include <cstddef>
#include <limits>

#include "PolyTree.hpp"

namespace pyscan {

	Point<> Segment::get_line() const {

		return Point<>(, this[1], this[2]);
	}



	bool Segment<>::lte(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(l);
            return alte(line[2], this->operator[](2) * a);
        } else {
            return line.above_closed(l_end_pt) && line.above_closed(r_end_pt);
        }
	}

    bool Segment<>::lt(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(l);
            return alt(line[2], this->operator[](2) * a);
        } else {
            return line.above(l_end_pt) && line.above(r_end_pt);
        }
    }


    bool Segment<>::gt(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(l);
            return alt(this->operator[](2) * a, line[2]);
        } else {
            return line.below(l_end_pt) && line.below(r_end_pt);
        }
    }

    bool Segment<>::gte(Point<> const& line) const {
        if (parallel(*this, line)) {
            double a = this->normalization_factor(l);
            return alte(this->operator[](2) * a, line[2]);
        } else {
            return line.below_closed(l_end_pt) && line.below_closed(r_end_pt);
        }
    }

    std::tuple<std::shared_ptr<PolyNode>, std::shared_ptr<PolyNode>>
    PolyNode::split(Point<> const& line) {
	    /*
	     * Line the left pt is the
	     */
	    Segment seg = line_to_segment(seg);
    	auto p1 = std::make_shared<PolyNode>();
    	auto p2 = std::make_shared<PolyNode>();

    	/*
    	 * Split the boundary segments
    	 */
        segment_list upper_boundary_seg;
        segment_list lower_boundary_seg;

        for (auto bs_it = boundary_seg.begin(); bs_it != boundary_seg.end(); ++bs_it) {
            if (seg.crossed(*bs_it)) {
                Segment upper, lower;
                std::tie(upper, lower) = bs_it->split(seg);
                lower_boundary_seg.push_back(lower);
                upper_boundary_seg.push_back(upper);
            } else if (seg.gte(*bs_it)) {
                upper_boundary_seg.push_back(*bs_it);
            } else {
                lower_boundary_seg.push_back(*bs_it);
            }
        }
        upper_boundary_seg.push_back(seg);
        lower_boundary_seg.push_back(seg);

        /*
         * split the crossing segments
         */
        wsegment_list upper_crossing_seg;
        wsegment_list lower_crossing_seg;

        for (auto bs_it = crossing_segments.begin(); bs_it != crossing_segments.end(); ++bs_it) {
            if (seg.crossed(*bs_it)) {
                WSegment upper, lower;
                std::tie(upper, lower) = bs_it->split(seg);
                lower_crossing_seg.push_back(lower);
                upper_crossing_seg.push_back(upper);
            } else if (seg.gte(*bs_it)) {
                upper_crossing_seg.push_back(*bs_it);
            } else {
                lower_crossing_seg.push_back(*bs_it);
            }
        }

        /*
         * Split the points.
         */

        line_list upper_pts;
        line_list lower_pts;
        for (auto p_it = points.begin(); p_it != points.end(); ++p_it) {
            if (seg.above_closed(*p_it)) {
                upper_pts.push_back(*p_it);
            } else {
                lower_pts.push_back(*p_it);
            }
        }

        upper_boundary_seg.push_back(seg);
        lower_boundary_seg.push_back(seg);
        return PolyNode(upper_boundary_seg, upper_weights, )
    }

    Segment line_to_segment(Point<2> const& ) const;
    Segment good_line_split() const;
    Segment good_vertex_split() const;

    ProjSegment PolyNode::good_line_split(int k) const {
        return weighted_choice(crossing_segments, weights);
    }

    ProjSegment PolyNode::good_vertex_split(int k) const;
    double PolyNode::get_weight() const;
}