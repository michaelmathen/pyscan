//
// Created by mmath on 2/12/18.
//

#ifndef PYSCAN_CUTTINGS_HPP
#define PYSCAN_CUTTINGS_HPP

#include <vector>
#include <memory>

namespace pyscan {

    class ProjObject {
        double a, b, c;

    public:
        ProjObject(double a, double b, double c) : a(a), b(b), c(c) {}

        ProjObject intersection(ProjObject const& other_line);

        bool below(ProjObject const& obj1);
        bool below_closed(ProjObject const& obj1);

        bool order(ProjObject const& l1, ProjObject const& l2);
        bool x_order(ProjObject const& p2);

    };


    class ProjSegment : public ProjObject {

    public:
        ProjSegment(ProjObject const& line, double xl, double xr) : ProjObject(line),
    };

    using node_ptr = std::shared_ptr<Node>;
    using line_list = std::vector<ProjObject>;
    using weight_list = std::vector<double>;

    class Node {


    };

    class SegmentNode : public Node {
        node_ptr upper;
        node_ptr lower;
        ProjObject line;

        //This is the x and scale parameter of the left and right
        // end of the segment.
        double x_l, s_l;
        double x_r, s_r;
    public:
        bool crossing(ProjObject const& line,)
    };

    class PolyNode : public Node {
        line_list boundary_lines;
        line_list crossing_lines;
        weight_list weights;
    public:
        PolyNode() {}

        PolyNode split(ProjObject line
    };


    std::vector<Cell> triangle_cutting(std::vector<ProjObject> const& lines, int r);


}
#endif //PYSCAN_CUTTINGS_HPP
