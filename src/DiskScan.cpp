//
// Created by mmath on 9/25/17.
//
#include <algorithm>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <assert.h>

#include "Statistics.hpp"
#include "Point.hpp"
#include "DiskScan.hpp"

namespace pyscan {

    template <typename W>
    void getLoc(Point<W, 2> pt, double& x1, double& x2) {
        x1 = get<0>(pt);
        x2 = get<1>(pt);
    }
    template <typename W>
    double getMeasured(Point<W, 2> const& pt) {
        return pt.getRedWeight();
    }

    template <typename W>
    double getBaseline(Point<W, 2> const& pt) {
        return pt.getBlueWeight();
    }


    template <typename W>
    double getX(Point<W, 2> const& pt) {
        return get<0>(pt);
    }


    template <typename W>
    double getY(Point<W, 2> const& pt) {
        return get<1>(pt);
    }

    template <typename W>
    bool colinear(Point<W> const& pt1,
                  Point<W> const& pt2,
                  Point<W> const& pt3){
        double x1, x2, x3, y1, y2, y3;
        getLoc(pt1, x1, y1);
        getLoc(pt2, x2, y2);
        getLoc(pt3, x3, y3);

        double
                a11 = x2 - x1, a12 = y2 - y1,
                a21 = x2 - x3, a22 = y2 - y3;
        return (a11 * a22 - a12 * a21 == 0);
    }

    template <typename W>
    void solveCircle3(Point<W> const& pt1, Point<W> const& pt2, Point<W> const& pt3,
                      double &a, double &b) {
        double x1, x2, x3, y1, y2, y3;
        getLoc(pt1, x1, y1);
        getLoc(pt2, x2, y2);
        getLoc(pt3, x3, y3);
        // Setup a matrix equation of the form Ax = b

        // A
        double a11 = x2 - x1, a12 = y2 - y1, a21 = x2 - x3, a22 = y2 - y3;
        // b
        double b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2, b2 = (y2 * y2
                                                                       + x2 * x2 - y3 * y3 - x3 * x3) / 2;

        double detA = a11 * a22 - a12 * a21;
        // inverse of A
        double ai11 = a22 / detA, ai12 = -a12 / detA, ai21 = -a21 / detA, ai22 = a11
                                                                                 / detA;

        // A^{-1} b = x
        a = ai11 * b1 + ai12 * b2;
        b = ai21 * b1 + ai22 * b2;
    }

    template <typename W>
    void solveCircle3(Point<W, 2> const& pt1, Point<W, 2> const& pt2, Point<W, 2> const& pt3,
                     double &a, double &b, double &r) {
        double x1, x2, x3, y1, y2, y3;
        getLoc(pt1, x1, y1);
        getLoc(pt2, x2, y2);
        getLoc(pt3, x3, y3);
        // Setup a matrix equation of the form Ax = b

        // A
        double a11 = x2 - x1, a12 = y2 - y1, a21 = x2 - x3, a22 = y2 - y3;
        // b
        double b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2, b2 = (y2 * y2
                                                                       + x2 * x2 - y3 * y3 - x3 * x3) / 2;

        double detA = a11 * a22 - a12 * a21;
        // inverse of A
        double ai11 = a22 / detA, ai12 = -a12 / detA, ai21 = -a21 / detA, ai22 = a11
                                                                                 / detA;

        // A^{-1} b = x
        a = ai11 * b1 + ai12 * b2;
        b = ai21 * b1 + ai22 * b2;
        r = sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));
    }


    template<typename W>
    bool sameLoc(Point<W> const& p1, Point<W> const& p2) {
        return getX(p1) == getX(p2) && getY(p1) == getY(p2);
    }

    bool onLineSegment(Point<double> const& pt1,
                       Point<double> const& pt2,
                       Point<double> const& pt3) {
        if (colinear(pt1, pt2, pt3)) {
            if (sameLoc(pt1, pt2) || sameLoc(pt2, pt3))
                return true;
            //Now we know that the point is on the same line
            if (getX(pt1) != getX(pt2)) {
                double theta = (getY(pt1) - getX(pt3)) / (getX(pt1) - getX(pt2));
                return (theta <= 1) && (theta >= 0);
            } else {
                double theta = (getY(pt1) - getY(pt3)) / (getY(pt1) - getY(pt2));
                return (theta <= 1) && (theta >= 0);
            }
        } else {
            return false;
        }
    }

    void findPerpVect(Point<double> const& p1, Point<double> const& p2, double* u, double* v) {
        double x1, x2, y1, y2;
        getLoc(p1, x1, y1);
        getLoc(p2, x2, y2);
        *u = y2 - y1;
        *v = x1 - x2;
    }

    inline void dot(double A[4], double B[2], double X[2]){
        X[0] = A[0] * B[0] + A[1] * B[1];
        X[1] = A[2] * B[0] + A[3] * B[1];
    }


    template<typename T, typename F, typename G>
    inline void partial_counts(T begin, T end,
                               std::vector<double> const& partitions,
                               std::vector<double>& counts, F& orderF, G& valueF) {
        //Partitions based on the break points.
        for (; begin != end; begin++) {
            auto lb = lower_bound(partitions.begin(), partitions.end(),
                                  orderF(*begin));
            counts[lb - partitions.begin()] += valueF(*begin);
        }
    }

    using val_t = std::tuple<size_t, double>;
    using crescent_t = std::vector<val_t>;

    template<typename T, typename F, typename G>
    inline void partial_counts_label(T begin, T end,
                               std::vector<double> const& partitions,
                               std::vector<crescent_t >& counts, F& orderF, G& valueF) {
        //Partitions based on the break points.
        for (; begin != end; begin++) {
            auto lb = lower_bound(partitions.begin(), partitions.end(),
                                  orderF(*begin));
            counts[lb - partitions.begin()].emplace_back(begin->getLabel(), valueF(*begin));
        }
    }


    double updateCounts(std::unordered_map<size_t, size_t>& curr_counts,
                        crescent_t& adding, crescent_t& removing) {
        double update_diff = 0;
        for (auto& val_pair : removing) {
            auto curr_el = curr_counts.find(std::get<0>(val_pair));
            //Start here.
            if (curr_el != curr_counts.end()) {
                if (curr_el->second <= 1) {
					update_diff -= std::get<0>(val_pair); 
                    curr_counts.erase(curr_el);
                } else {
                    curr_counts[curr_el->first] = curr_el->second + 1;
                }
            } else {
                assert("The current set contains element that are being removed in the set");
            }
        }
        for (auto& val_pair : adding) {
            auto curr_el = curr_counts.find(std::get<0>(val_pair));
            //Start here.
            if (curr_el != curr_counts.end()) {
                curr_counts[curr_el->first] = curr_el->second + 1;
            } else {
				update_diff += std::get<0>(val_pair); 
                curr_counts[std::get<0>(val_pair)] = 1;
            }
        }
		return update_diff;
    }

    template< typename T, typename F>
    double computeLabelTotal(T begin, T end, F func, std::unordered_map<size_t, size_t>& label_map) {
        double total = 0;
        for (; begin != end; ++begin) {
            if (label_map.end() == label_map.find(begin->getLabel())) {
                total += func(*begin);
                label_map[begin->getLabel()] = 0;
            }
            label_map[begin->getLabel()] = label_map[begin->getLabel()] + 1;
        }
        return total;
    }

    template< typename T, typename F>
    double computeLabelTotal(T begin, T end, F func) {
        std::unordered_map<size_t, size_t> label_map;
        return computeLabelTotal(begin, end, func, label_map);
    }

    template <typename F>
    Disk diskScanLabels(std::vector<LPoint<double, 2>>& net,
                        std::vector<LPoint<double, 2>>& sampleM,
                        std::vector<LPoint<double, 2>>& sampleB, F scan) {
        /*
         * Each point is defined with a label. For each division we store the set of
         * labels.
         */
        //Calculate the total measured and baseline value.
        auto msBegin = sampleM.begin();
        auto bsBegin = sampleB.begin();
        auto msEnd = sampleM.end();
        auto bsEnd = sampleB.end();

        auto nB = net.begin();
        auto nE = net.end();
        std::vector <Point<double>> netSampleSorted(nB, nE);

        Disk currMax;
        double maxStat = 0;
        //Need to change this so that it keeps track of the number of counts we have
        //seen.


        std::vector<crescent_t> mCountsA(nE - nB, crescent_t());
        std::vector<crescent_t> bCountsA(nE - nB, crescent_t());
        std::vector<crescent_t> mCountsR(nE - nB, crescent_t());
        std::vector<crescent_t> bCountsR(nE - nB, crescent_t());
        double m_Total = computeLabelTotal(msBegin, msEnd, getMeasured<double>);
        double b_Total = computeLabelTotal(bsBegin, bsEnd, getBaseline<double>);

        auto sortedB = netSampleSorted.begin();
        auto sortedEOuter = netSampleSorted.end();
        for (auto i = nB; i != nE - 2; i++) {
            sortedEOuter = std::remove(sortedB, sortedEOuter, *i);
            auto sortedE = sortedEOuter;
            for (auto j = i + 1; j != nE - 1; j++) {
                auto el = std::find(sortedB, sortedE, *j);
                std::swap(*el, *(sortedE - 1));
                sortedE = sortedE - 1;
                //Create a vector between the two points
                double orthoX, orthoY;
                findPerpVect(*i, *j, &orthoX, &orthoY);
                double cX = (get<0>(*i) + get<0>(*j)) / 2.0;
                double cY = (get<1>(*i) + get<1>(*j)) / 2.0;
                auto isNotCol = [&i, &j](Point<double> const& pt) {
                    double x1, x2, x3, y1, y2, y3;
                    getLoc(*i, x1, y1);
                    getLoc(*j, x2, y2);
                    getLoc(pt, x3, y3);
                    double
                            a11 = x2 - x1, a12 = y2 - y1,
                            a21 = x2 - x3, a22 = y2 - y3;
                    return (a11 * a22 - a12 * a21 != 0);
                };


                // Partition these into a set of adding points and removing points
                auto partitionF = [orthoX, orthoY, cX, cY](Point<double> const &pt) {
                    return (get<0>(pt) - cX) * orthoX + (get<1>(pt) - cY) * orthoY <= 0;
                };

                auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Point<double> const &pt) {
                    // If the point lines up with either of the reference
                    // point then we take this to be a disk defined by only
                    // the reference points.
                    // We are projecting a vector created between
                    //the disk center and center point between the two points.
                    double a, b;
                    solveCircle3(*i, *j, pt, a, b);
                    return orthoX * (a - cX) + orthoY * (b - cY);
                };
                auto compF = [&orderF](Point<double> const &pt1, Point<double> const &pt2) {
                    return orderF(pt1) < orderF(pt2);
                };

                auto asIterEnd = std::partition(msBegin, msEnd, isNotCol);
                auto bsIterEnd = std::partition(bsBegin, bsEnd, isNotCol);
                auto nIterEnd = std::partition(sortedB, sortedE, isNotCol);

                double mCount = 0;
                double bCount = 0;
                // will have two sets an added set and a current set.
                // added set is for stuff that will never leave. The current set on the other hand
                // will correspond to every trajectory that we overlap.

                std::unordered_map<size_t, size_t> m_curr_set;
                std::unordered_map<size_t, size_t> b_curr_set;

                std::for_each(asIterEnd, msEnd, [&](LPoint<double> const &pt) {
                    if (onLineSegment(*i, *j, pt)) {
                        if (m_curr_set.find(pt.getLabel()) == m_curr_set.end()) {
                            mCount += getMeasured(pt);
                            m_curr_set[pt.getLabel()] = 1;
                        } else {
                            m_curr_set[pt.getLabel()] = m_curr_set[pt.getLabel()] + 1;
                        }
                    }
                });
                std::for_each(bsIterEnd, bsEnd, [&](LPoint<double> const &pt) {
                    if (onLineSegment(*i, *j, pt)) {
                        if (b_curr_set.find(pt.getLabel()) == b_curr_set.end()) {
                            bCount += getBaseline(pt);
                            b_curr_set[pt.getLabel()] = 1;
                        } else {
                            b_curr_set[pt.getLabel()] = b_curr_set[pt.getLabel()] + 1;
                        }
                    }
                });

                std::sort(sortedB, nIterEnd, compF);
                std::vector<double> orderV;
                orderV.reserve(nIterEnd - sortedB);
                for (auto b = sortedB; b != nIterEnd; b++) {
                    orderV.push_back(orderF(*b));
                }
                auto mHigherIt = std::partition(msBegin, asIterEnd, partitionF);
                auto bHigherIt = std::partition(bsBegin, bsIterEnd, partitionF);

                std::fill(mCountsR.begin(), mCountsR.end(), crescent_t());
                std::fill(bCountsR.begin(), bCountsR.end(), crescent_t());
                std::fill(mCountsA.begin(), mCountsA.end(), crescent_t());
                std::fill(bCountsA.begin(), bCountsA.end(), crescent_t());
                /*Probably most of the time is spent here*/
                partial_counts_label(msBegin, mHigherIt, orderV, mCountsR, orderF, getMeasured<double>);
                mCount = mCount + computeLabelTotal(msBegin, mHigherIt, getMeasured<double>, m_curr_set);
                partial_counts_label(bsBegin, bHigherIt, orderV, bCountsR, orderF, getBaseline<double>);
                bCount = bCount + computeLabelTotal(bsBegin, bHigherIt, getBaseline<double>, b_curr_set);
                partial_counts_label(mHigherIt, asIterEnd, orderV, mCountsA, orderF,
                               getMeasured<double>);
                partial_counts_label(bHigherIt, bsIterEnd, orderV, bCountsA, orderF,
                               getBaseline<double>);
                /*----------------------------------------------*/
                //Now scan over the counts.
                auto size = nIterEnd - sortedB;
                for (int k = 0; k < size; k++) {
					mCount += updateCounts(m_curr_set, mCountsA[k], mCountsR[k]);
					bCount += updateCounts(b_curr_set, bCountsA[k], bCountsR[k]);
                    double m_hat = mCount / m_Total;
                    double b_hat = bCount / b_Total;
                    double newStat = scan(m_hat, b_hat);
                    if (maxStat <= newStat) {
                        double a, b, r;
                        solveCircle3(*i, *j, *(sortedB + k), a, b, r);
                        Disk currDisk(newStat, a, b, r);
                        currMax = currDisk;
                        maxStat = newStat;
                    }
                }
            }
        }
        return currMax;
    }

    template <typename F>
    Disk diskScan(std::vector<Point<double, 2>>& net, std::vector<Point<double, 2>>& sampleM, std::vector<Point<double, 2>>& sampleB, F scan) {

        //Calculate the total measured and baseline value.
        double m_Total = 0;
        double b_Total = 0;
        auto msBegin = sampleM.begin();
        auto bsBegin = sampleB.begin();
        auto msEnd = sampleM.end();
        auto bsEnd = sampleB.end();
        for (auto aIt = msBegin; aIt != msEnd; aIt++)
            m_Total += getMeasured(*aIt);
        for (auto bIt = bsBegin; bIt != bsEnd; bIt++)
            b_Total += getBaseline(*bIt);

        auto nB = net.begin();
        auto nE = net.end();
        std::vector <Point<double>> netSampleSorted(nB, nE);

        Disk currMax;
        double maxStat = 0;
        std::vector<double> mCountsA(nE - nB, 0);
        std::vector<double> bCountsA(nE - nB, 0);
        std::vector<double> mCountsR(nE - nB, 0);
        std::vector<double> bCountsR(nE - nB, 0);

        auto sortedB = netSampleSorted.begin();
        auto sortedEOuter = netSampleSorted.end();
        for (auto i = nB; i != nE - 2; i++) {
            sortedEOuter = std::remove(sortedB, sortedEOuter, *i);
            auto sortedE = sortedEOuter;
            for (auto j = i + 1; j != nE - 1; j++) {
                auto el = std::find(sortedB, sortedE, *j);
                std::swap(*el, *(sortedE - 1));
                sortedE = sortedE - 1;
                //Create a vector between the two points
                double orthoX, orthoY;
                findPerpVect(*i, *j, &orthoX, &orthoY);
                double cX = (get<0>(*i) + get<0>(*j)) / 2.0;
                double cY = (get<1>(*i) + get<1>(*j)) / 2.0;
                auto isNotCol = [&i, &j](Point<double> const& pt) {
                    double x1, x2, x3, y1, y2, y3;
                    getLoc(*i, x1, y1);
                    getLoc(*j, x2, y2);
                    getLoc(pt, x3, y3);
                    double
                        a11 = x2 - x1, a12 = y2 - y1,
                        a21 = x2 - x3, a22 = y2 - y3;
                    return (a11 * a22 - a12 * a21 != 0);
                };



          // Partition these into a set of adding points and removing points
          auto partitionF = [orthoX, orthoY, cX, cY](Point<double> const &pt) {
              return (get<0>(pt) - cX) * orthoX + (get<1>(pt) - cY) * orthoY <= 0;
          };

          auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Point<double> const &pt) {
              // If the point lines up with either of the reference
              // point then we take this to be a disk defined by only
              // the reference points.
              // We are projecting a vector created between
              //the disk center and center point between the two points.
              double a, b;
              solveCircle3(*i, *j, pt, a, b);
              return orthoX * (a - cX) + orthoY * (b - cY);
          };
          auto compF = [&orderF](Point<double> const &pt1, Point<double> const &pt2) {
              return orderF(pt1) < orderF(pt2);
          };

          auto asIterEnd = std::partition(msBegin, msEnd, isNotCol);
          auto bsIterEnd = std::partition(bsBegin, bsEnd, isNotCol);
          auto nIterEnd = std::partition(sortedB, sortedE, isNotCol);

          double mCount = 0;
          double bCount = 0;
          std::for_each(asIterEnd, msEnd, [i, j, &mCount](Point<double> const &pt) {
              if (onLineSegment(*i, *j, pt)) {
                mCount += getMeasured(pt);
              }
          });
          std::for_each(bsIterEnd, bsEnd, [i, j, &bCount](Point<double> const &pt) {
              if (onLineSegment(*i, *j, pt)) {
                bCount += getMeasured(pt);
              }
          });

          std::sort(sortedB, nIterEnd, compF);
          std::vector<double> orderV;
          orderV.reserve(nIterEnd - sortedB);
          for (auto b = sortedB; b != nIterEnd; b++) {
            orderV.push_back(orderF(*b));
          }
          auto mHigherIt = std::partition(msBegin, asIterEnd, partitionF);
          auto bHigherIt = std::partition(bsBegin, bsIterEnd, partitionF);

          std::fill(mCountsR.begin(), mCountsR.end(), 0);
          std::fill(bCountsR.begin(), bCountsR.end(), 0);
          std::fill(mCountsA.begin(), mCountsA.end(), 0);
          std::fill(bCountsA.begin(), bCountsA.end(), 0);
          /*Probably most of the time is spent here*/
          partial_counts(msBegin, mHigherIt, orderV, mCountsR, orderF, getMeasured<double>);
          mCount = mHigherIt - msBegin + mCount;
          partial_counts(bsBegin, bHigherIt, orderV, bCountsR, orderF, getBaseline<double>);
          bCount = bHigherIt - bsBegin + bCount;
          partial_counts(mHigherIt, asIterEnd, orderV, mCountsA, orderF,
                         getMeasured<double>);
          partial_counts(bHigherIt, bsIterEnd, orderV, bCountsA, orderF,
                         getBaseline<double>);
          /*----------------------------------------------*/
          //Now scan over the counts.
          auto size = nIterEnd - sortedB;
          for (int k = 0; k < size; k++) {
            mCount += mCountsA[k] - mCountsR[k];
            bCount += bCountsA[k] - bCountsR[k];
            double m_hat = mCount / m_Total;
            double b_hat = bCount / b_Total;
            double newStat = scan(m_hat, b_hat);
            if (maxStat <= newStat) {
                double a, b, r;
                solveCircle3(*i, *j, *(sortedB + k), a, b, r);
                Disk currDisk(newStat, a, b, r);
                currMax = currDisk;
                maxStat = newStat;
            }
          }
        }
      }
      return currMax;
    }


    Disk diskScanStatLabels(std::vector<LPoint<double, 2>>& net, std::vector<LPoint<double, 2>>& sampleM, std::vector<LPoint<double, 2>>& sampleB, double rho) {
        return diskScanLabels(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }

    Disk diskScanStat(std::vector<Point<double, 2>>& net, std::vector<Point<double, 2>>& sampleM, std::vector<Point<double, 2>>& sampleB, double rho) {
        return diskScan(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }
}
