//
// Created by mmath on 2/12/18.
//
#include <tuple>
#include <algorithm>
#include <cstddef>
#include <limits>

#include "Cuttings.hpp"

namespace pyscan {

	double det2(double a1, double b1, double a2, double b2) {
		return a1 * b2 - b1 * a2;
	}


	bool approx_zero(double val) {
		return fabs(val) <= std::numeric_limits<double>::epsilon();
	}

	bool approx_lte_zero(double val) {
		return val <= std::numeric_limits<double>::epsilon();
	}

	bool approx_lt_zero(double val) {
		return val < -std::numeric_limits<double>::epsilon();
	}

	bool approx_gte_zero(double val) {
		return val >= std::numeric_limits<double>::epsilon();
	}

	ProjObject ProjObject::intersection(ProjObject const& other_line) const {
		return ProjObject(det2(b, c, other_line.b, other_line.c), 
							det2(a, c, other_line.a, other_line.c),
							det2(a, b, other_line.a, other_line.b)
			);
	}

    bool ProjObject::dual_is_below_closed(ProjObject const& pt) const {
    	return approx_lte_zero(this->evaluate(pt));
    }
  	
  	bool ProjObject::dual_is_above_closed(ProjObject const& pt) const {
  		return approx_gte_zero(this->evaluate(pt));
  	}

  	double ProjObject::evaluate(ProjObject const& pt) {
  		return a * pt.a + b * pt.b + c * pt.c;
  	}

  	bool ProjObject::parallel(ProjObject const& line) const {
  		return approx_zero(det2(a, b, line.a, line.b));
  	}
	
	ProjObject ProjSegment::get_line() const {
		return ProjObject(a, b, c);
	}

    bool ProjSegment::above_closed(ProjObject const& other_line) const {
    	if (this->parallel(other_line)) {
            return approx_gte_zero(this->c - other_line->c);
    	} 
    	return other_line.dual_is_below_closed(l_end_point) && 
    			other_line.dual_is_below_closed(r_end_point);
    }

    bool ProjSegment::below_closed(ProjObject const& other_line) const {
        if (this->parallel(other_line)) {
            return approx_lte_zero(this->c - other_line->c);
        } 
        return other_line.dual_is_above_closed(l_end_point) && 
                other_line.dual_is_above_closed(r_end_point);
    }

    bool ProjSegment::is_crossed(ProjObject const& line) const {
        return (other_line.dual_is_above(l_end_point) && 
                other_line.dual_is_below(r_end_point)) || 
               (other_line.dual_is_below(l_end_point) && 
                other_line.dual_is_above(r_end_point));
    }

 `

    bool SegmentNode::below_closed(ProjSegment const& segment) const {
    	return segment.below_closed(this->line);
    }

    bool SegmentNode::above_closed_pt(ProjObject const& pt) const {
    	return line.dual_is_above_closed(pt);
    }

    bool SegmentNode::below_closed_pt(ProjObject const& pt) const {
    	return line.dual_is_below_closed(pt);
    }

    void split_lines(segment_list const& segs, weight_list const& weights, ProjObject const& lines, 
        segment_list& upper, weight_list& u_ws, segment_list& lower, weight_list& l_ws) {

    }

    void split_points(line_list const& segs, ProjObject const& lines, line_list& upper, line_list& lower) {

    }

    void PolyNode::split(ProjObject const& line, node_ptr& p1, node_ptr& p2) const {
    	p1 = std::make_shared<PolyNode>();
    	p2 = std::make_shared<PolyNode>();

        segment_list upper_boundary_segs;
        segment_list lower_boundary_segs;
        split_lines(this->boundary_segments, line, upper_boundary_segs, lower_boundary_segs);

        segment_list upper_segs, lower_segs;
        weight_list u_w, l_w;
        split_lines(this->crossing_segments, line, upper_segs, u_w, lower_segs, l_w);

        line_list upper_pts, lower_pts;
        split_points(this->points, line, upper_pts, lower_pts);
    }

    ProjSegment PolyNode::good_line_split(int k) const {
        return weighted_choice(crossing_segments, weights);
    }

    ProjSegment PolyNode::good_vertex_split(int k) const;
    double PolyNode::get_weight() const;
}