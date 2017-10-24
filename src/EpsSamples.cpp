//
// Created by mmath on 7/26/17.
//

#include <limits>
#include <deque>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <random>
#include <cassert>
#include <boost/functional/hash.hpp>

#include "EpsSamples.hpp"


namespace pyscan {


    double yVal(Line2 l1, double xval) {
        return std::get<0>(l1) * xval + std::get<1>(l1);
    }

    double yVal(Line l1, double xval) {
        return std::get<0>(l1) * xval + std::get<1>(l1);
    }

    template<typename... T>
    struct Hash {
        size_t operator()(std::tuple<T...> const& arg) const noexcept {
            return boost::hash_value(arg);
        }
    };

    double xLoc(double top_a, double top_b, double a2, double b2) {
        // b1 or b2 could be infinite
        return (b2 - top_b) / (top_a - a2);
    }

    double xLoc(Line2 l1, Line2 l2) {
        return xLoc(std::get<0>(l1), std::get<1>(l1), std::get<0>(l2), std::get<1>(l2));
    }


    Line2 loc(Line2 l1, Line2 l2) {
        double x = xLoc(std::get<0>(l1), std::get<1>(l1), std::get<0>(l2), std::get<1>(l2));
        double y = yVal(l1, x);
        return Line2(x, y);
    }
//
//   bool Edge::crosses(Line2 const &crossing_line) const {
//       double crossing_pt = xLoc(std::get<0>(crossing_line), std::get<1>(crossing_line),
//            std::get<0>(edge_line), std::get<1>(edge_line));
//
//       return x_begin <= crossing_pt && crossing_pt <= x_end;
//   }
//
//   bool Edge::sameLine(Edge const& other_edge) const {
//       return std::get<0>(edge_line) == std::get<0>(other_edge.edge_line) &&
//               std::get<1>(edge_line) == std::get<1>(other_edge.edge_line);
//   }
//
//   void split_edge(edge_ptr edge, Line2 const& intersect_line, edge_ptr e1, edge_ptr e2) {
//       // Split the second crossing edge
//       double crossing_pt = xLoc(std::get<0>(intersect_line), std::get<1>(intersect_line),
//                                 std::get<0>(edge->edge_line), std::get<1>(edge->edge_line));
//
//       e1 = std::make_shared<Edge>(*edge, edge->x_begin, crossing_pt);
//       e2 = std::make_shared<Edge>(*edge, crossing_pt, edge->x_end);
//   }
//
//   bool Vertex::crosses(Line2 const &line) const {
//       return std::get<0>(line) * std::get<0>(vertex_location) + std::get<1>(line) == std::get<1>(vertex_location);
//   }
//
//
//   Vertex::Vertex(Line2 const &l1, Line2 const &l2, simplex_ptr pl, simplex_ptr pr) : left_ptr(pl), right_ptr(pr) {
//       double x_value = xLoc(std::get<0>(l1), std::get<1>(l1), std::get<0>(l2), std::get<1>(l2));
//       vertex_location = Line2(x_value, x_value * std::get<0>(l1) + std::get<1>(l1));
//   }
//   /*
//   Simplex Simplex::getNext(Line2 const& new_line) {
//
//   }
//   */
//   Line2 toLine2(Line const& l) {
//       return Line2(std::get<0>(l), std::get<1>(l));
//   }
//
//   //////////////////////////////////////////////////////
//   //Splitters///////////////////////////////////////////
//   //////////////////////////////////////////////////////
//
//   double Splitter::bottom() const {
//       if (_unbounded && !_standing) {
//           return -std::numeric_limits<double>::infinity();
//       } else {
//           double x_val = getX();
//           return std::min(yVal(l, x_val), yVal(termination, x_val));
//       }
//   }
//
//   double Splitter::top() const {
//       if (_unbounded && _standing) {
//           return std::numeric_limits<double>::infinity();
//       } else {
//           double x_val = getX();
//           return std::max(yVal(l, x_val), yVal(termination, x_val));
//       }
//   }
//
//   double Splitter::getX() const {
//       return xLoc(l, r);
//   }
//
//   bool Splitter::crosses(Line2 const& nl) const {
//       double x_val = getX();
//       double line_y_val = yVal(nl, x_val);
//       return bottom() <= line_y_val && line_y_val <= top();
//   }
//
//   bool Splitter::crosses(Line2 const& nl, Line2& pt) const {
//       /*
//        * Computes the crossing point and stores it in the pt;
//        */
//       double x_val = getX();
//       double line_y_val = yVal(nl, x_val);
//       pt = Line2(x_val, line_y_val);
//       return bottom() <= line_y_val && line_y_val <= top();
//
//   }
//
//   bool Splitter::crosses(Line const& nl) const {
//       return crosses(toLine2(nl));
//   }
//
//   bool Splitter::passesBelow(Line2 const& nl) const {
//       //Checks to see if the line passes below the bottom terminating vertex
//       if (_unbounded && !_standing) {
//           return false;
//       } else {
//           double x_val = getX();
//           return yVal(nl, x_val) <= bottom();
//       }
//   }
//   bool Splitter::passesAbove(Line2 const& nl) const {
//       //Checks to see if the line passes above the bottom terminating vertex
//       if (_unbounded && _standing) {
//           return false;
//       } else {
//           double x_val = getX();
//           return top() <= yVal(nl, x_val);
//       }
//   }
//
//   bool Splitter::leftOf(Line2 const& nl) const {
//       return rightOf(nl);
//   }
//
//   bool Splitter::rightOf(Line2 const& nl) const {
//       //If it is unbounded then it has to go through
//       // either l or r line
//       if (passesAbove(nl) != _standing) {
//           //Check to see if it crosses the left line to the left
//           // of the splitter
//           return xLoc(nl, l) <= getX();
//       } else {
//           //It must cross the term line either to the left or right
//           //of the splitter
//           return xLoc(nl, termination) <= getX();
//       }
//   }
//
//   bool Splitter::leftOf(Line const& nl) const {
//       return leftOf(toLine2(nl));
//   }
//   bool Splitter::rightOf(Line const& nl) const {
//       return rightOf(toLine2(nl));
//   }
//
//   Line2 Splitter::crossPt() const {
//       return Line2(getX(), yVal(l, getX()));
//   }
//
//   bool Splitter::isCrossLine(Line2 const& nl) const {
//       return nl == l || nl == r;
//   }
//
//
//   /*
//    *
//    * Copies the entire top half into the new simplex. If the simplex is open then entrance or exit could be -1 then
//    * this also preserves the location of the opening. This does not adjust the edge pointers.
//    */
//   void Simplex::copyTopHalf(int entrance, int exit,
//                             edge_ptr entrance_edge,
//                             edge_ptr exit_edge,
//                             edge_ptr cutting_edge,
//                             vertex_ptr entrance_v,
//                             vertex_ptr exit_v,
//                             std::vector<edge_ptr>& top_edges,
//                             std::vector<vertex_ptr>& top_vertices) {
//       if (entrance == -1) {
//           top_edges.push_back(cutting_edge);
//           top_edges.push_back(exit_edge);
//           top_edges.insert(top_edges.end(), edge_list.begin() + exit + 1, edge_list.end());
//           top_vertices.push_back(exit_v);
//           top_vertices.insert(top_vertices.end(), corner_ptrs.begin() + exit, corner_ptrs.end());
//       } else if (exit == -1) {
//           top_edges.insert(top_edges.end(), edge_list.begin(), edge_list.begin() + entrance);
//           top_edges.push_back(entrance_edge);
//           top_edges.push_back(cutting_edge);
//           top_vertices.insert(top_vertices.end(), corner_ptrs.begin(), corner_ptrs.begin() + entrance);
//           top_vertices.push_back(entrance_v);
//       } else if (entrance < exit) {
//           top_edges.insert(top_edges.end(), edge_list.begin(), edge_list.begin() + entrance);
//           top_edges.push_back(entrance_edge);
//           top_edges.push_back(cutting_edge);
//           top_edges.push_back(exit_edge);
//           top_edges.insert(top_edges.end(), edge_list.begin() + exit + 1, edge_list.end());
//
//           top_vertices.insert(top_vertices.end(), corner_ptrs.begin(), corner_ptrs.end()+ entrance);
//           top_vertices.push_back(entrance_v);
//           top_vertices.push_back(exit_v);
//           top_vertices.insert(top_vertices.end(), top_vertices.begin() + exit, top_vertices.end());
//
//       } else {
//           top_edges.push_back(exit_edge);
//           top_edges.insert(top_edges.end(), edge_list.begin() + exit + 1, edge_list.begin() + entrance);
//           top_edges.push_back(entrance_edge);
//           top_edges.push_back(cutting_edge);
//
//           top_vertices.insert(top_vertices.end(), corner_ptrs.begin() + exit + 1, corner_ptrs.end() + entrance);
//           top_vertices.push_back(entrance_v);
//           top_vertices.push_back(exit_v);
//       }
//   }
//
//   void Simplex::copyBottomHalf(int entrance, int exit,
//                                edge_ptr entrance_edge,
//                                edge_ptr exit_edge,
//                                edge_ptr cutting_edge,
//                                vertex_ptr entrance_v,
//                                vertex_ptr exit_v,
//                                std::vector<edge_ptr>& top_edges,
//                                std::vector<vertex_ptr>& top_vertices) {
//      copyTopHalf(exit, entrance, exit_edge, entrance_edge, cutting_edge, exit_v, entrance_v, top_edges, top_vertices);
//   }
//
//   void Simplex::updatePtrs(simplex_ptr old_ptr) {
//       /*
//        * For each pointer it either points inside or outside of the cell. Since we
//        * created a new cell in an earlier step we need to adjust the pointer to point
//        * at ourselves
//        */
//       for (auto& edge : edge_list) {
//           if (edge->top_ptr.lock() == old_ptr) {
//               edge->top_ptr = this_ptr;
//           } else {
//               edge->bottom_ptr = this_ptr;
//           }
//       }
//       for (auto& vertex : corner_ptrs) {
//           if (vertex->right_ptr.lock() == old_ptr) {
//               vertex->right_ptr = this_ptr;
//           } else {
//               vertex->left_ptr = this_ptr;
//           }
//       }
//   }
//
//   bool above(Line l, Line2 vertex) {
//       return std::get<0>(l) * std::get<0>(vertex) + std::get<1>(l) > std::get<1>(vertex);
//   }
//
//   bool below(Line l, Line2 vertex) {
//       return std::get<0>(l) * std::get<0>(vertex) + std::get<1>(l) < std::get<1>(vertex);
//   }
//
//   void seperateLines(LineList const& lines, Line2 left_vertex, Line2 right_vertex,
//                      LineList above_list, LineList below_list){
//       for (Line l : lines) {
//           //This should put the line in the top or bottom set.
//           //It might not end up in either if it is the same line as the seperator.
//           if (above(l, left_vertex) || above(l, right_vertex)) {
//               above_list.push_back(l);
//           }
//           if (below(l, left_vertex) || below(l, right_vertex)) {
//               below_list.push_back(l);
//           }
//       }
//   }
//
//   void Simplex::computeLRVertexCross(simplex_ptr left_ptr, Line2 l, vertex_ptr left_vt, vertex_ptr right_vt) {
//       /*
//        * Check the line and see if it crosses any vertices. The left entrance and right entrance could potentialy
//        * cross vertices. If it crosses a vertex then we will set left_vt to the left vertex or right_vt to the right
//        * vertex. If the left or right entrance is not crossed then we set them to nullptrs.
//        */
//       vertex_ptr first_ptr = nullptr, second_ptr = nullptr;
//       for (auto vertex : corner_ptrs) {
//           if (vertex->crosses(l) && first_ptr == nullptr) {
//               first_ptr == vertex;
//           } else if (vertex->crosses(l)) {
//               second_ptr = vertex;
//           }
//       }
//
//       if (first_ptr == nullptr || left_ptr == nullptr || left_ptr != first_ptr->left_ptr.lock()) {
//           left_vt = second_ptr, right_vt = first_ptr;
//       } else {
//           left_vt = first_ptr, right_vt = second_ptr;
//       }
//   }
//
//   void Simplex::computeLRcrossing(simplex_ptr left_ptr, Line2 l, edge_ptr left_edge, edge_ptr right_edge) const {
//       /*
//        * Computes the left and right edge crossings. If there is no left or no right crossing then we set
//        * left or right to false.
//        */
//       edge_ptr first_cross = nullptr, second_cross = nullptr;
//       for (auto e : this->edge_list) {
//           if (e->crosses(l) && first_cross == nullptr) {
//               if (first_cross == nullptr)
//                   first_cross = e;
//               else
//                   second_cross = e;
//           }
//       }
//       // Three cases fc, sc null, sc null, neither null.
//       if (first_cross == nullptr && second_cross == nullptr) {
//           left_edge = nullptr, right_edge = nullptr;
//       } else if (second_cross == nullptr) {
//           if (left_ptr != nullptr) {
//               left_edge = first_cross;
//           } else {
//               right_edge = first_cross;
//           }
//       } else {
//           if (first_cross->bottom_ptr.lock() == left_ptr ||
//                   first_cross->top_ptr.lock() == left_ptr) {
//               left_edge = first_cross;
//               right_edge = second_cross;
//           } else {
//               left_edge = second_cross;
//               right_edge = first_cross;
//           }
//       }
//   }
//
//   edge_ptr Simplex::findOppositeEdge(edge_ptr intersecting_edge, double x_val) const {
//       for (edge_ptr & e : this->edge_list) {
//           if (e->x_begin <= x_val && e->x_end <= x_val && e != intersecting_edge) {
//               return e;
//           }
//       }
//       return nullptr;
//   }
//
//   void insertSplitter(Splitter const& spl, std::vector<LineList> & crossing_list, std::vector<Splitter>& splitter_list) {
//       /*
//        * Inserts the splitter into the splitter list and modifies the crossing_list by splitting all lines in
//        * the crossing list into the set that lies on the left and right sides of the crossing list.
//        */
//       auto s_it = std::find_if(splitter_list.begin(), splitter_list.end(), [&](Splitter const& sp){
//           return spl.getX() < sp.getX();
//       });
//       splitter_list.insert(s_it, spl);
//       LineList left_lines;
//       LineList right_lines;
//       for (auto line : crossing_list[s_ix]){
//           if (spl.crosses(line)) {
//               left_lines.push_back(line);
//               right_lines.push_back(line);
//           } else if (spl.leftOf(line)) { //splitter is to the left of the line
//               right_lines.push_back(line);
//           } else {
//               left_lines.push_back(line);
//           }
//       }
//       //Insert the new line crossings into the crossing and splitter list
//       auto s_ix = s_it - splitter_list.begin();
//       crossing_list[s_ix] = right_lines;
//       crossing_list.insert(crossing_list.begin() + s_ix, left_lines);
//
//   }
//
//   Splitter Simplex::makeSplitter(edge_ptr e1, Line2 new_line, edge_ptr term) {
//       if (term == nullptr) {
//           return {e1->edge_line, new_line, e1->top_ptr.lock() == this};
//       } else {
//           return {e1->edge_line, new_line, term->edge_line};
//       }
//   }
//
//   void mergeSplitters(Line2 const& new_line, bool above, std::vector<Splitter>& splitters, std::vector<LineList>& crossings) {
//       /*
//        * Enumerate the splitters and merge cells that should not belong.
//        */
//       std::set<Line> cell_elements;
//       std::vector<Splitter> new_splitters;
//       std::vector<LineList> new_crossings;
//
//       for (size_t i = 0; i < splitters.size(); i++) {
//           //Need to handle the other cases.
//           if (splitters[i].crosses(new_line)) {
//               cell_elements.insert(crossings[i].begin(), crossings[i].end());
//
//               if (above != splitters[i].standing()) {
//                   new_crossings.emplace_back(cell_elements.begin(), cell_elements.end());
//                   new_splitters.emplace_back(splitters[i], new_line);
//                   cell_elements = std::set<Line>();
//               }
//           } else if (splitters[i].passesAbove(new_line) == above) {
//               new_crossings.emplace_back(crossings[i].begin(), crossings.end());
//               new_splitters.push_back(splitters[i]);
//           }
//       }
//       // insert the last cell.
//       cell_elements.insert(crossings[splitters.size() + 1].begin(), crossings[splitters.size() + 1].end());
//       new_crossings.emplace_back(cell_elements.begin(), cell_elements.end());
//       splitters = new_splitters;
//       crossings = new_crossings;
//   }
//
//   void splitSimplex(simplex_ptr& sp, simplex_ptr& entr_ptr, simplex_ptr& upper, simplex_ptr& lower){
//
//   }
//
//   void Simplex::splitSplitters(Line2 const& new_line, simplex_ptr entr_ptr,
//                                decltype(crossings)& above_crossings, decltype(splitter_list)& above_splitters,
//                                decltype(corner_ptrs)& above_corners,
//                                decltype(crossings)& below_crossings, decltype(splitter_list)& below_splitters,
//                                decltype(corner_ptrs)& below_corners) {
//       /*
//        * Takes splitters and removes lines that fall no longer fall between two splitters.
//        * Since they earlier fell inside, but now we only need to check whether the line
//        * passes above the two end points of adjacent splitters.
//        */
//       // Always will be at least one region.
//       size_t i = 0;
//       edge_ptr left_edge, right_edge;
//       vertex_ptr  left_vt, right_vt;
//       computeLRVertexCross(std::move(entr_ptr), new_line, left_vt, right_vt);
//       computeLRcrossing(std::move(entr_ptr), new_line, left_edge, right_edge);
//       decltype(crossings) local_crossings = crossings;
//       decltype(splitter_list) local_splitters = splitter_list;
//       decltype(corner_ptrs) local_corners = corner_ptrs;
//       // (1) Insert the new splitters and local crossings into the list and
//       // also add the newly created edges and (possibly vertices).
//       // TODO handle the case where the line goes through a vertex.
//       if (left_vt == nullptr && left_edge != nullptr) {
//           auto term_left = findOppositeEdge(left_edge, xLoc(left_edge->edge_line, new_line));
//           //In these cases the simplex is open on the top or the bottom so we need to extend
//           //either down or up. The simplex is convex so we find out which side it opens to.
//           Splitter left_splitter = makeSplitter(left_edge, new_line, term_left);
//           insertSplitter(left_splitter, local_crossings, local_splitters);
//       }
//       if (right_vt == nullptr && right_edge != nullptr) {
//           auto term_right = findOppositeEdge(right_edge, xLoc(right_edge->edge_line, new_line));
//           Splitter right_splitter = makeSplitter(right_edge, new_line, term_right);
//           insertSplitter(right_splitter, local_crossings, local_splitters);
//       }
//
//       // (2) Now need to cut all regions defined by splitters into an upper and lower portion depending on
//       // how the new line cuts them.
//       // If the line passes below the region then crossings and splitters end up in the above
//       // If the line pass above the region then crossings and splitters end up in the below
//       // If the line passes through the region then crossings and splitters end up in both.
//       if (splitter_list.size() == 0) {
//           above_crossings.push_back(local_crossings[0]);
//           below_crossings.push_back(local_crossings[0]);
//       } else {
//           // We have at least 1 splitter and 2 crossings
//           // We check each splitter in sequence until we start getting crossings.
//           // We will cross the first splitter and last splitters we inserted
//           Line2 left_entrance, right_entrance;
//           size_t i = 0;
//           for (; i < local_splitters.size(); i++) {
//               if (local_splitters[i].passesBelow(new_line)) {
//                   above_crossings.push_back(local_crossings[i]);
//                   above_splitters.push_back(local_splitters[i]);
//               } else if (local_splitters[i].passesAbove(new_line) {
//                   below_crossings.push_back(local_crossings[i]);
//                   below_splitters.push_back(local_splitters[i]);
//               }
//               if (local_splitters[i].isCrossLine(new_line)){
//                   //Now have left entrance
//                   left_entrance = local_splitters[i].crossPt();
//                   i++; //Set to next splitter
//                   break;
//               }
//           }
//           // Now everything between.
//           LineList mergeIntoBelow;
//           LineList mergeIntoAbove;
//           for (; i < local_splitters.size(); i++) {
//               if (local_splitters[i].crosses(new_line, right_entrance)){
//
//                   LineList above_lines, below_lines;
//                   seperateLines(local_crossings[i], left_entrance, right_entrance, above_lines, below_lines);
//                   above_crossings.push_back(above_lines);
//                   below_crossings.push_back(below_lines);
//                   above_splitters.push_back(local_splitters[i]);
//                   below_splitters.push_back(local_splitters[i]);
//                   /*
//                   if (local_splitters[i].standing()) {
//                       below_splitters.emplace_back(Splitter(local_splitters[i], new_line));
//                       mergeIntoAbove.insert(mergeIntoAbove.end(), local_crossings[i].begin(), local_crossings[i].end());
//                       below_crossings.push_back()
//                   } else {
//                       above_splitters.emplace_back(Splitter(local_splitters[i], new_line));
//                       mergeIntoBelow.insert(mergeIntoAbove.end(), local_crossings[i].begin(), local_crossings[i].end());
//                   }
//                   */
//                   //Now have left entrance
//                   left_entrance = right_entrance;
//               } else if (local_splitters[i].passesBelow(new_line)) {
//                   above_crossings.push_back(local_crossings[i]);
//                   above_splitters.push_back(local_splitters[i]);
//               } else if (local_splitters[i].passesAbove(new_line) {
//                   below_crossings.push_back(local_crossings[i]);
//                   below_splitters.push_back(local_splitters[i]);
//               }
//           }
//       }
//       // Now need to merge splitters.
//       mergeSplitters(new_line, true, above_splitters, above_crossings);
//       mergeSplitters(new_line, false, below_splitters, below_crossings);
//   }
//
//   Simplex Simplex::upper_split(Line2 const& new_line, edge_ptr clock_wise, vertex_ptr lower_simplex) {
//       /*
//        * Computes the upper split of a simplex
//        */
//
//       Simplex top_simplex;
//
//       // 4 cases for closed simplex, v -> e, e -> e, e -> v, and v -> v
//       // 4 cases for open simplex, v -> o, o -> v, o-> e, e -> o.
//
//       bool entrance_vertex;
//       bool exit_vertex;
//       if (entrance_vertex && exit_vertex) {
//
//       } else if (entrance_vertex && !exit_vertex) {
//
//       } else if (!entrance_vertex && exit_vertex) {
//
//       } else {
//
//       }
//       return Simplex();
//   }
//
//   bool Simplex::closed() const {
//       return !open;
//   }
//
//
//   Simplex Simplex::split(Line2 const& new_line,
//                          edge_ptr const& clock_entrance_edge,
//                          edge_ptr const& counter_entrance_edge,
//                          vertex_ptr upper_simplex,
//                          vertex_ptr lower_simplex
//   ) {
//
//       return Simplex();
//   }

   void splitCell(Cell const &cell, Line const &line, std::deque<Cell> &cells, std::vector<Cell> &output);

    Trapezoid::Trapezoid(double ta, double tb, double lx, double ba, double bb, double rx) : top_a(ta),
                                                                                             top_b(tb),
                                                                                             bottom_a(ba),
                                                                                             bottom_b(bb),
                                                                                             left_x(lx),
                                                                                             right_x(rx) {}
    Trapezoid::Trapezoid() : top_a(0),
                             top_b(std::numeric_limits<double>::infinity()),
                             bottom_a(0),
                             bottom_b(-std::numeric_limits<double>::infinity()),
                             left_x(-std::numeric_limits<double>::infinity()),
                             right_x(std::numeric_limits<double>::infinity()) {}

/*
 * This does not consider incident lines to be crossing lines.
 */
    bool Trapezoid::crossing(double a, double b) const {
        /*
         * Checks to see if this line crosses a corner or
         * side.
         */
        return (crossesTop(a, b) + crossesBottom(a, b) + crossesRightSide(a, b) + crossesLeftSide(a, b)) > 1;
        /*
        if ((a == top_a && b == top_b) || (a == bottom_a && b == bottom_b)) {
            return false; // line is the same as the top or bottom line.
        } else {
            double lefty = a * left_x + b,
                    righty = a * right_x + b;
            double lty = top_a * left_x + top_b,
                    rty = top_a * right_x + top_b,
                    lby = bottom_a * left_x + bottom_b,
                    rby = bottom_a * right_x + bottom_b;
            if (lefty >= lty) { // crosses above the left most top point, but crosses bellow the right most point.
                return righty <= rty;
            } else if (lefty >= lby) { //goes through left side.
                return !(lefty == lby && righty < rby); // Not (If it goes through the corner
                //but passes below the rightmost point)
            } else {
                // In this case it goes below the left bottom point so just need to check if it crosses
                // above the rightmost bottom point.
                return righty >= rby;
            }
        }
        */

    }

    bool Trapezoid::incident(double a, double b) const {
        /*
         * Checks to see if this just crosses the boundary of the trapezoid
         */
        if ((top_a == a && top_b == b) || (bottom_a == a && b == bottom_b)) {
            return true;
        } else {

            double top_ly = left_x * top_a + top_b;
            double top_ry = right_x * top_a + top_b;
            double bottom_ly = left_x * bottom_a + bottom_b;
            double bottom_ry = right_x * bottom_a + bottom_b;
            double ly = left_x * a + b;
            double ry = right_x * a + b;
            // Checks to see if the one corner is the same, but the other passes above.
            // If bottom_ry or top_ry have infinite magnitude then the strictly greater than
            // should prevent incidence.
            return (ly == top_ly && ry > top_ry) ||
                    (ly == bottom_ly && ry < bottom_ry) ||
                    (ry == bottom_ry && ly < bottom_ly) ||
                    (ry == top_ry && ly > top_ly);
        }
    }


    void Trapezoid::setBottomLine(double a, double b) {
        /*
        if (bottom_a * left_x + bottom_b > a * left_x + b) {
            std::cout << "Problem happened: Bad Bottom" << std::endl;
        }
         */
        this->bottom_b = b;
        this->bottom_a = a;
    }

    void Trapezoid::setTopLine(double a, double b) {
        /*
        if (bottom_a * left_x + bottom_b > top_b * left_x + top_b) {
            std::cout << "Problem happened: Bad Top" << std::endl;
        }
         */
        this->top_a = a;
        this->top_b = b;
    }

    void Trapezoid::setLeft(double left) {
        if (left == std::numeric_limits<double>::infinity()) {
            std::cout << "Set Left is going bad" << std::endl;
        }
        this->left_x = left;
    }

    void Trapezoid::setRight(double right) {
        if (right == -std::numeric_limits<double>::infinity()) {
            std::cout << "Set Right is going bad" << std::endl;
        }
        this->right_x = right;
    }

    /*
    * These consider incident lines to be crossing lines
    */
    bool Trapezoid::crossesTop(double a, double b) const {
        double lty = left_x * top_a + top_b;
        double rty = top_a * right_x + top_b;
        double ly = a * left_x + b;
        double ry = a * right_x + b;
        if (right_x == std::numeric_limits<double>::infinity() &&
            left_x == -std::numeric_limits<double>::infinity()) {
            return true;
        } else if (right_x == std::numeric_limits<double>::infinity()) {
            if (ly == lty) return true;
            else if (ly > lty && a < top_a) return true;
            else if (ly < lty && a > top_a) return true;
            else return false;
        } else if (left_x == -std::numeric_limits<double>::infinity()) {
            if (ry == rty) return true;
            else if (ry > rty && a > top_a) return true;
            else if (ry < rty && a < top_a) return true;
            else return false;
        } else if (top_b == std::numeric_limits<double>::infinity()){
            return false; // finite width interval, but infinitely tall trapezoid
        } else {
            return (ly <= lty && rty <= ry) != (ly >= lty && rty >= ry);
        }
    }

    bool Trapezoid::crossesBottom(double a, double b) const {
        double lby = left_x * bottom_a + bottom_b;
        double rby = bottom_a * right_x + bottom_b;
        double ly = a * left_x + b;
        double ry = a * right_x + b;
        if (right_x == std::numeric_limits<double>::infinity() &&
            left_x == -std::numeric_limits<double>::infinity()) {
            return true;
        } else if (right_x == std::numeric_limits<double>::infinity()) {
            if (ly == lby) return true;
            else if (ly > lby && a < bottom_a) return true;
            else if (ly < lby && a > bottom_a) return true;
            else return false;
        } else if (left_x == -std::numeric_limits<double>::infinity()) {
            if (ry == rby) return true;
            else if (ry > rby && a > bottom_a) return true;
            else if (ry < rby && a < bottom_a) return true;
            else return false;
        } else if (bottom_b == -std::numeric_limits<double>::infinity()){
            return false; // finite width interval, but infinitely tall trapezoid
        } else {
            return (ly <= lby && rby <= ry) != (ly >= lby && rby >= ry);
        }
    }

    bool Trapezoid::crossesLeftSide(double a, double b) const {
        if (top_b == std::numeric_limits<double>::infinity() &&
                bottom_b == -std::numeric_limits<double>::infinity()) {
            return true;
        } else if (left_x == -std::numeric_limits<double>::infinity()) {
            if (top_b == std::numeric_limits<double>::infinity()) {
                return a < bottom_a;
            } else if (bottom_b == -std::numeric_limits<double>::infinity()) {
                return a > top_a;
            } else {
                return a < bottom_a && a > top_a;
            }
        } else  {
            double lty = top_a * left_x + top_b;
            double lby = bottom_a * left_x + bottom_b;
            double ly = a * left_x + b;
            return lty >= ly && ly >= lby;
        }
    }

    bool Trapezoid::crossesRightSide(double a, double b) const {
        if (top_b == std::numeric_limits<double>::infinity() &&
            bottom_b == -std::numeric_limits<double>::infinity()) {
            return true;
        } else if (right_x == std::numeric_limits<double>::infinity()) {
            if (top_b == std::numeric_limits<double>::infinity()) {
                return a > bottom_a;
            } else if (bottom_b == -std::numeric_limits<double>::infinity()) {
                return a < top_a;
            } else {
                return a > bottom_a && a < top_a;
            }
        } else  {
            double rty = top_a * right_x + top_b;
            double rby = bottom_a * right_x + bottom_b;
            double ry = a * right_x + b;
            return rty >= ry && ry >= rby;
        }
    }

    bool Trapezoid::finiteTop() const {
        return top_b != std::numeric_limits<double>::infinity();
    }
    bool Trapezoid::finiteBottom() const {
        return bottom_b != -std::numeric_limits<double>::infinity();
    }
    bool Trapezoid::finiteLeft() const {
        return left_x != -std::numeric_limits<double>::infinity();
    }
    bool Trapezoid::finiteRight() const {
        return right_x != std::numeric_limits<double>::infinity();
    }

    bool Trapezoid::isDegenerate() const {
        return this->degenerate != TrapType::NORMAL;
    }
    Line2 Trapezoid::urCorner() const {
        return Line2(right_x, top_a * right_x + top_b);
    }
    Line2 Trapezoid::llCorner() const {
        return Line2(left_x, bottom_a * left_x + bottom_b);
    }
    Line2 Trapezoid::lrCorner() const {
        return Line2(right_x, bottom_a * right_x + bottom_b);
    }
    Line2 Trapezoid::ulCorner() const {
        return Line2(left_x, top_a * left_x + top_b);
    }
    bool Trapezoid::lDegenerate() const {
        return this->degenerate == TrapType::L_Deg;
    }

    bool Trapezoid::rDegenerate() const {
        return this->degenerate == TrapType::R_Deg;
    }

    std::string Trapezoid::print() const {
        std::stringstream os;
        os << "Trapezoid(" << top_a << ", " << top_b << ", " << left_x << ", "
           << bottom_a << ", " << bottom_b << ", " << right_x << ")";
        return os.str();
    }

    bool Trapezoid::contains(Line2 const& pt) const {
        return std::get<0>(pt) * top_a + top_b >= std::get<1>(pt) &&
               left_x <= std::get<0>(pt) &&
               std::get<0>(pt) <= right_x &&
               std::get<1>(pt) >= std::get<0>(pt) * bottom_a + bottom_b;
    }

    Cell computeCrossings(Trapezoid const &trap, Cell const& cell) {
        long weight = 0;
        LineList new_line_list;
        PointList new_point_list;
        for (auto line : std::get<1>(cell)) {
            if (trap.crossesInterior(std::get<0>(line), std::get<1>(line))) {
                new_line_list.push_back(line);
                weight += 1 << std::get<2>(line);
            }
        }
        for (auto point : std::get<2>(cell)) {
            if (trap.contains(point)) {
                new_point_list.push_back(point);
            }
        }
        return Cell(trap, new_line_list, new_point_list, weight);
    }


    void splitCell(Cell const &cell, Line const &line, std::vector<Cell> &cells) {
        /*
         * Takes a single trapezoid and splits it into 2 to 4 trapezoids depending on the
         * orientation. Determines what lines cross each of the resulting trapezoids and
         * then adds them to the cells or output depending on their size.
         *
         * This is not handling degenerate cases right now...
         */
        auto &trap = std::get<0>(cell);
        double a = std::get<0>(line),
                b = std::get<1>(line);
        if (trap.crossesLeftSide(a, b) && trap.crossesRightSide(a, b)) {
            /*
             *  Could also cross the top and bottom if this is an infinite sized trapezoid
             */

            Trapezoid new_trap1 = trap, new_trap2 = trap;
            LineList crossing_lines1, crossing_lines2;
            new_trap1.bottom_a = a, new_trap1.bottom_b = b;
            new_trap2.top_a = a, new_trap2.top_b = b;

            cells.emplace_back(computeCrossings(new_trap1, cell));
            cells.emplace_back(computeCrossings(new_trap2, cell));

        } else if (trap.crossesRightSide(a, b) != trap.crossesLeftSide(a, b)) {
            /*
             * Both the top and the bottom need to be finite in this case.
             */
            Trapezoid new_trap1 = trap,
                    new_trap2 = trap,
                    new_trap3 = trap;
            double x_ins_pt;
            if (trap.crossesTop(a, b) && trap.finiteTop()) {
                x_ins_pt = xLoc(trap.top_a, trap.top_b, a, b);
                if (trap.crossesRightSide(a, b)) {
                    new_trap3.setLeft(x_ins_pt);
                    new_trap2.setTopLine(a, b);
                } else {
                    new_trap3.setRight(x_ins_pt);
                    new_trap1.setTopLine(a, b);
                }
                new_trap3.setBottomLine(a, b);
            } else {
                x_ins_pt = xLoc(trap.bottom_a, trap.bottom_b, a, b);
                if (trap.crossesRightSide(a, b)) {
                    new_trap3.setLeft(x_ins_pt);
                    new_trap2.setBottomLine(a, b);
                } else {
                    new_trap3.setRight(x_ins_pt);
                    new_trap1.setBottomLine(a, b);
                }
                new_trap3.setTopLine(a, b);
            }
            new_trap1.setRight(x_ins_pt);
            new_trap2.setLeft(x_ins_pt);
            new_trap3.degenerate = trap.crossesRightSide(a, b) ? TrapType::L_Deg : TrapType::R_Deg;
            LineList crossing_lines1, crossing_lines2, crossing_lines3;
            cells.emplace_back(computeCrossings(new_trap1, cell));
            cells.emplace_back(computeCrossings(new_trap2, cell));
            cells.emplace_back(computeCrossings(new_trap3, cell));
        } else if (trap.crossesTop(a, b) && trap.crossesBottom(a, b)) {


            Trapezoid new_trap1 = trap,
                    new_trap2 = trap,
                    new_trap3 = trap,
                    new_trap4 = trap;
            double rx = std::max(xLoc(trap.top_a, trap.top_b, a, b), xLoc(trap.bottom_a, trap.bottom_b, a, b));
            double lx = std::min(xLoc(trap.top_a, trap.top_b, a, b), xLoc(trap.bottom_a, trap.bottom_b, a, b));
            new_trap1.setRight(lx);
            new_trap2.setLeft(lx);
            new_trap2.setRight(rx);
            new_trap3.setLeft(lx);
            new_trap3.setRight(rx);
            new_trap4.setLeft(rx);
            if (a > 0) {
                new_trap2.setTopLine(a, b);
                new_trap3.setBottomLine(a, b);
                new_trap2.degenerate = TrapType::R_Deg;
                new_trap3.degenerate = TrapType::L_Deg;
            } else {
                new_trap3.setBottomLine(a, b);
                new_trap2.setTopLine(a, b);
                new_trap2.degenerate = TrapType::L_Deg;
                new_trap3.degenerate = TrapType::R_Deg;
            }
            cells.emplace_back(computeCrossings(new_trap1, cell));
            cells.emplace_back(computeCrossings(new_trap2, cell));
            cells.emplace_back(computeCrossings(new_trap3, cell));
            cells.emplace_back(computeCrossings(new_trap4, cell));
        } else {
            std::cout << "got here somehow " << std::get<0>(cell).print() << std::endl;
            std::cout << std::get<0>(line) << " " << std::get<1>(line) << std::endl;
        }
    }

    bool Trapezoid::crossesInterior(double a, double b) const {
        return this->crossing(a, b) && !this->incident(a, b);
    }
    bool randomGreedyCut(LineList &lines, PointList &points, int r, Cell & output) {

        /*
         * Start with all lines in a single infinite sized trapezoid.
         * Randomly add lines and split any trapezoid that is crossed by too many lines.
         * Proceed till all trapezoids are crossed by less than r lines.
         * Returns a set of trapezoids with a list of lines that cross each one.
         * This should take approximately O(nr) time.
         */
        if (r == 1) {
            return false;
        }
        long weight = 0;
        std::for_each(lines.begin(), lines.end(), [&](Line& el){
            weight += 1 << std::get<2>(el);
        });
        std::vector<Cell> cells;
        cells.emplace_back(Cell(Trapezoid(), lines, points, weight));
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(lines.begin(), lines.end(), gen);
        for (auto &line : lines) {
            std::vector<Cell> new_cells;
            for (auto &cell : cells) {
                if (std::get<3>(cell) <= weight / r) {
                    //std::cout << std::get<3>(cell) << " " << std::get<2>(cell).size() << std::endl;
                    if (std::get<2>(cell).size() >= points.size() / (r * r)) {
                        output = cell;
                        return true;
                    }
                } else if (std::get<0>(cell).crossesInterior(std::get<0>(line), std::get<1>(line))) {
                    splitCell(cell, line, new_cells);
                } else {
                    new_cells.push_back(cell);
                }

            }
            cells = new_cells;
            if (cells.empty()) {
                break;
            }
        }
        return false;
    }

    std::vector<Cell> randomCuttings(LineList & lines, PointList &points, int r) {
        /*
         * Start with all lines in a single infinite sized trapezoid.
         * Randomly add lines and split any trapezoid that is crossed by too many lines.
         * Proceed till all trapezoids are crossed by less than r lines.
         * Returns a set of trapezoids with a list of lines that cross each one.
         * This should take approximately O(nr) time.
         */
        if (r == 1) {
            return std::vector<Cell>{Cell(Trapezoid(), lines, points, 0)};
        }
        long weight = 0;
        std::for_each(lines.begin(), lines.end(), [&](Line& el){
            weight += 1 << std::get<2>(el);
        });
        std::vector<Cell> cells;
        cells.emplace_back(Cell(Trapezoid(), lines, points, weight));
        std::vector<Cell> non_active_cells;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(lines.begin(), lines.end(), gen);
        for (auto &line : lines) {
            std::vector<Cell> new_cells;
            for (auto &cell : cells) {
                if (std::get<3>(cell) <= weight / r) {
                    non_active_cells.push_back(cell);
                } else if (std::get<0>(cell).crossesInterior(std::get<0>(line), std::get<1>(line))) {
                    splitCell(cell, line, new_cells);
                } else {
                    new_cells.push_back(cell);
                }

            }
            cells = new_cells;
            if (cells.empty()) {
                break;
            }
        }
        return non_active_cells;
    }

    Line toDual(Line const& line) {
        return Line(std::get<0>(line), -std::get<1>(line), std::get<2>(line));
    }
    Line2 toDual(Line2 const& line) {
        return Line2(std::get<0>(line), -std::get<1>(line));
    }

    LineList dualizeToLines(std::vector<Cell> const& cells) {
        /*
         * Prunes out intersections that do not belong to the set of lines.
         */
        std::unordered_set<Line2, Hash<double, double> > unique_lines;
        auto addWrapper = [&](Line2 const& line) {
            if (std::get<0>(line) == -std::numeric_limits<double>::infinity() ||
                    std::get<0>(line) == std::numeric_limits<double>::infinity() ||
                    std::get<1>(line) == -std::numeric_limits<double>::infinity() ||
                    std::get<1>(line) == std::numeric_limits<double>::infinity()) {
                return;
            }
            unique_lines.emplace(line);
        };
        for (auto cell : cells) {
            auto trap = std::get<0>(cell);
            addWrapper(toDual(Line2(trap.ulCorner())));
            addWrapper(toDual(Line2(trap.urCorner())));
            addWrapper(toDual(Line2(trap.llCorner())));
            addWrapper(toDual(Line2(trap.lrCorner())));
        }
        LineList lines;
        for (auto line : unique_lines) {
            lines.emplace_back(Line(std::get<0>(line), std::get<1>(line), 1));
        }

        return lines;
    }

    LineList dualizeToLines(std::vector<Line2> const& points) {
        LineList lines;
        for(auto line : points) {
            lines.emplace_back(Line(std::get<0>(line), std::get<1>(line), 1));
        }
        return lines;
    }

    Line doubleWeight(Line const& line) {
        return Line(std::get<0>(line), std::get<1>(line), std::get<2>(line) + 1);
    }
    /*
    void scaleDownCell(Cell& old_cell, double s) {

        auto cell = std::get<0>(old_cell);
        auto pts = std::get<2>(old_cell);
        auto pt_value = [&] (Line2 const& p1) {
            return std::get<0>(p1) * cell.top_a + cell.top_b - std::get<1>(p1);
        };
        auto pt_order = [&](Line2 const& l1, Line2 const& l2) {
            return pt_value(l1) > pt_value(l2);
        };

        std::(pts.begin(), pts.end(), pt_order);


    }
    */
    /*
    Cell createCell(PointList const& points, LineList const& lines) {
        // Create trapezoid containing these points.
        auto min_max_pair = std::minmax_element(points.begin(), points.end(), [](Line2 p1, Line2 p2){
            return std::get<0>(p1) < std::get<0>(p2);
        });

    }
    */

    std::vector<Cell> createPartition(PointList & points, int r) {

        /*
         * Generate a test set for this set of points.
         */
        double s = points.size() / (r * r); // size of cells.
        PointList empty_points;
        auto lines = dualizeToLines(points);
        auto cells = randomCuttings(lines, empty_points, r);
        auto test_set = dualizeToLines(cells);

        std::vector<Cell> partitions;

        for (int i = 0; i < r * r; i++) {

            /*
             * Cut the test set into trapezoids based on the current weighting.
             */
            if (points.size() <= points.size() / (r * r) ) {
                //partitions.emplace_back(createCell(points, test_set));
                return partitions;
            }
            int r_i = static_cast<int>(round(sqrt(points.size() / s)) + .5);
            //std::cout << "r_i " << r_i << " points " << points.size() << std::endl;
            Cell cell;
            if (!randomGreedyCut(test_set, points, r_i, cell)) {
                //std::cout << "couldn't find a partition" << std::endl;
                return partitions;
            }
            // Move one of the sides down until the region contains exactly s points.
            //scaleDownCell(cell, s);
            //std::cout << "found a partition"  << std::endl;
            partitions.push_back(cell);

            /*
             * Remove the points from the set.
             */
            std::unordered_set<Line2, Hash<double, double> > points_set(std::get<2>(cell).begin(), std::get<2>(cell).end());
            PointList new_point_set;
            for (Line2 pt : points) {
                if (points_set.find(pt) == points_set.end()) {
                    new_point_set.push_back(pt);
                }
            }
            points = new_point_set;

            /*
             * Double the weight of any lines that cross the cell found in the triangulation.
             */
            std::unordered_set<Line, Hash<double, double, int> > line_set(std::get<1>(cell).begin(), std::get<1>(cell).end());
            LineList updated_test_set;
            for (Line line : test_set) {
                if (line_set.find(line) == line_set.end()) {
                    updated_test_set.push_back(line);
                } else {
                    updated_test_set.emplace_back(doubleWeight(line));
                }
            }
            test_set = updated_test_set;
        }
        return partitions;
    }


}

/*
void pyscan::halfpsaceSample(pyscan::point_it begin, pyscan::point_it end,
                             pyscan::point_it output_b, pyscan::point_it output_e) {

}
 */
