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
    		// (-by - c, ay, a)
    		// (b x, -a x - c, b)

    		// (-c, 0, a)
    		// (0, -c, b)
    		ProjObject test_pt;
    		if (approx_zero(a)) {
    			test_pt = ProjObject(0, -c, b);
    		} else {
    			test_pt = ProjObject(-c, 0, a);
    		}
    		return other_line.dual_is_below_closed(test_pt);
    	} 
    	return other_line.dual_is_below_closed(l_end_point) && 
    			other_line.dual_is_below_closed(r_end_point);
    }

    bool ProjSegment::below_closed(ProjObject const& other_line) const {
    	if (this->parallel(other_line)) {
    		// (-by - c, ay, a)
    		// (b x, -a x - c, b)

    		// (-c, 0, a)
    		// (0, -c, b)
    		ProjObject test_pt;
    		if (approx_zero(a)) {
    			test_pt = ProjObject(0, -c, b);
    		} else {
    			test_pt = ProjObject(-c, 0, a);
    		}
    		return other_line.dual_is_above_closed(test_pt);
    	} 
    	return other_line.dual_is_above_closed(l_end_point) && 
    			other_line.dual_is_above_closed(r_end_point);
    }

    bool ProjSegment::is_crossed(ProjObject const& line) const {
    	double l_v = line.evaluate(l_end_point);
    	double r_v = line.evaluate(r_end_point);
    	if (approx_zero(l_v)  || approx_zero(r_v)) {
    		return false;
    	} else {
    		return l_v * r_v < 0;
    	}
    }

  	bool SegmentNode::crossing(ProjSegment const& segment) const {
  		return segment.is_crossed(this->line);
  	}

  	bool SegmentNode::above_closed(ProjSegment const& segment) const {
  		return segment.above_closed(this->line);
  	}

    bool SegmentNode::below_closed(ProjSegment const& segment) const {
    	return segment.below_closed(this->line);
    }

    bool SegmentNode::above_closed_pt(ProjObject const& pt) const {
    	return line.dual_is_above_closed(pt);
    }
    bool SegmentNode::below_closed_pt(ProjObject const& pt) const {
    	return line.dual_is_below_closed(pt);
    }

    void PolyNode::split(ProjObject const& line, node_ptr& p1, node_ptr& p2) const {
    	p1 = std::make_shared<PolyNode>();
    	p2 = std::make_shared<PolyNode>();

    	line_list upper_boundary;
    	line_list lower_boundary;
    	for (auto& l : boundary_lines) {
    	}

    	for (auto& l : crossing_lines) {
    		
    	}
    }

    ProjSegment PolyNode::good_line_split(int k) const;
    ProjSegment PolyNode::good_vertex_split(int k) const;
    double PolyNode::get_weight() const;
}