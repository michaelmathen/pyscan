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

    void solveCircle3(Point<> const& pt1, Point<> const& pt2, Point<> const& pt3,
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

    void solveCircle3(Point<> const& pt1, Point<> const& pt2, Point<> const& pt3,
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

    void findPerpVect(Point<> const& p1, Point<> const& p2, double* u, double* v) {
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
            auto lb = std::lower_bound(partitions.begin(), partitions.end(),
                                  orderF(*begin));
            counts[lb - partitions.begin()] += valueF(*begin);
        }
    }


    template<typename T, typename F, typename G>
    inline void partial_counts_label(T begin, T end,
                               std::vector<double> const& partitions,
                               std::vector<crescent_t >& counts, F orderF, G valueF) {
        //Partitions based on the break points.
        for (; begin != end; begin++) {
            auto lb = std::lower_bound(partitions.begin(), partitions.end(),
                                  orderF(*begin));
            counts[lb - partitions.begin()].emplace_back(begin->getLabel(), valueF(*begin));
        }
    }

    double updateCounts(std::unordered_map<size_t, size_t>& curr_counts,
                        crescent_t& adding, crescent_t& removing) {
        /*
        * The adding and removing crescents define the
        */
        double update_diff = 0;
        for (auto& val_pair : adding) {
            auto curr_el = curr_counts.find(val_pair.label);
            //Start here.
            if (curr_el != curr_counts.end()) {
                curr_counts[curr_el->first] = curr_el->second + 1;
            } else {
                update_diff += val_pair.value;
                curr_counts[val_pair.label] = 1;
            }
        }
        for (auto& val_pair : removing) {
            auto curr_el = curr_counts.find(val_pair.label);
            //Start here.
            if (curr_el != curr_counts.end()) {
                if (curr_el->second <= 1) {
                    update_diff -= val_pair.value;
                    curr_counts.erase(curr_el);
                } else {
                    curr_counts[curr_el->first] = curr_el->second - 1;
                }
            } else {
                assert("The current set contains element that are being removed in the set");
            }
        }
        return update_diff;
    }

    template <typename F>
    Disk diskScanLabels(std::vector<LPoint<>>& net,
                        std::vector<LPoint<>>& sampleM,
                        std::vector<LPoint<>>& sampleB, F scan) {
        /*
         * Each point is defined with a label. For each division we store the set of
         * labels.
         */
        //Calculate the total measured and baseline value.
        auto msBegin = sampleM.begin();
        auto bsBegin = sampleB.begin();
        auto msEnd = sampleM.end();
        auto bsEnd = sampleB.end();
        double m_Total = computeLabelTotal(msBegin, msEnd, getMeasured);
        double b_Total = computeLabelTotal(bsBegin, bsEnd, getBaseline);

        auto nB = net.begin();
        auto nE = net.end();
        std::vector <LPoint<>> netSampleSorted(nB, nE);

        Disk currMax;
        double maxStat = 0;

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
                auto isNotCol = [&i, &j](LPoint<> const& pt) {
                    return !colinear(*i, *j, pt);
                };
                // Partition these into a set of adding points and removing points
                auto partitionF = [orthoX, orthoY, cX, cY](LPoint<> const &pt) {
                    return (get<0>(pt) - cX) * orthoX + (get<1>(pt) - cY) * orthoY <= 0;
                };

                auto orderF = [orthoX, orthoY, &i, &j, cX, cY](LPoint<> const &pt) {
                    // If the point lines up with either of the reference
                    // point then we take this to be a disk defined by only
                    // the reference points.
                    // We are projecting a vector created between
                    //the disk center and center point between the two points.
                    double a, b;
                    solveCircle3(*i, *j, pt, a, b);
                    return orthoX * (a - cX) + orthoY * (b - cY);
                };
                auto compF = [&orderF](LPoint<> const &pt1, LPoint<> const &pt2) {
                    return orderF(pt1) < orderF(pt2);
                };

                auto asIterEnd = std::partition(msBegin, msEnd, isNotCol);
                auto bsIterEnd = std::partition(bsBegin, bsEnd, isNotCol);
                auto nIterEnd = std::partition(sortedB, sortedE, isNotCol);

                // will have two sets an added set and a current set.
                // added set is for stuff that will never leave. The current set on the other hand
                // will correspond to every trajectory that we overlap.
                auto onSegment = [&](LPoint<> const& pt) {
                    return onLineSegment(*i, *j, pt);
                };

                //Counts the points that lie on the line segment between i and j.
                //These points are colinear so have been removed from the scan.
                std::unordered_map<size_t, size_t> m_curr_set;
                std::unordered_map<size_t, size_t> b_curr_set;
                double m_count = computeLabelTotalF(asIterEnd, msEnd, getMeasured, m_curr_set, onSegment);
                double b_count = computeLabelTotalF(bsIterEnd, bsEnd, getBaseline, b_curr_set, onSegment);

                std::sort(sortedB, nIterEnd, compF);
                std::vector<double> orderV;
                orderV.reserve(static_cast<size_t>(nIterEnd - sortedB));
                for (auto b = sortedB; b != nIterEnd; b++) {
                    orderV.push_back(orderF(*b));
                }

                //Partition the point set into an adding and removing sets.
                auto mHigherIt = std::partition(msBegin, asIterEnd, partitionF);
                auto bHigherIt = std::partition(bsBegin, bsIterEnd, partitionF);

                std::vector<crescent_t> mCountsA(static_cast<size_t>(nE - nB), crescent_t());
                std::vector<crescent_t> bCountsA(static_cast<size_t>(nE - nB), crescent_t());
                std::vector<crescent_t> mCountsR(static_cast<size_t>(nE - nB), crescent_t());
                std::vector<crescent_t> bCountsR(static_cast<size_t>(nE - nB), crescent_t());

                /*Probably most of the time is spent here*/
                partial_counts_label(msBegin, mHigherIt, orderV, mCountsR, orderF,
                                     getMeasured);
                m_count += computeLabelTotal(msBegin, mHigherIt, getMeasured, m_curr_set);
                partial_counts_label(bsBegin, bHigherIt, orderV, bCountsR, orderF,
                                     getBaseline);
                b_count += computeLabelTotal(bsBegin, bHigherIt, getBaseline, b_curr_set);
                partial_counts_label(mHigherIt, asIterEnd, orderV, mCountsA, orderF,
                               getMeasured);
                partial_counts_label(bHigherIt, bsIterEnd, orderV, bCountsA, orderF,
                               getBaseline);
                /*----------------------------------------------*/
                //Now scan over the counts.
                auto size = nIterEnd - sortedB;
                for (int k = 0; k < size; k++) {
                    m_count += updateCounts(m_curr_set, mCountsA[k], mCountsR[k]);
                    b_count += updateCounts(b_curr_set, bCountsA[k], bCountsR[k]);
                    double m_hat = m_count / m_Total;
                    double b_hat = b_count / b_Total;
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

    bool contains(double a, double b, double r, Point<> const& pt) {
        double la = (pt.getX() - a), lb = (pt.getY() - b);
        return la * la + lb * lb <= r * r;
    }

    template <typename F>
    Disk diskScanSlow(std::vector<Point<>>& net, std::vector<Point<>>& sampleM, std::vector<Point<>>& sampleB, F scan) {
        double m_Total = 0;
        double b_Total = 0;
        for (auto& aIt : sampleM)
            m_Total += getMeasured(aIt);
        for (auto& bIt : sampleB)
            b_Total += getBaseline(bIt);

        Disk maxDisk(0, 0, 0, 0);
        for (auto i = net.begin(); i != net.end() - 2; i++) {
            for (auto j = i + 1; j != net.end() - 1; j++) {
                for (auto k = j + 1; k != net.end(); k++) {
                    if (!colinear(*i, *j, *k)) {
                        double m_curr = 0;
                        double b_curr = 0;
                        double a, b, r;
                        solveCircle3(*i, *j, *k, a, b, r);
                        for (auto& p : sampleM) {
                           if (contains(a, b, r, p)) {
                              m_curr += getMeasured(p);
                           }
                        }
                        for (auto& p : sampleB) {
                            if (contains(a, b, r, p)) {
                                b_curr += getBaseline(p);
                            }
                        }
                        if (scan(m_curr / m_Total, b_curr / b_Total) >= maxDisk.fValue()) {
                            maxDisk = Disk(scan(m_curr / m_Total, b_curr / b_Total), a, b, r);
                        }

                    }
                }
            }
        }
        return maxDisk;
    }

    template <typename F>
    Disk diskScanSlowLabels(std::vector<LPoint<>>& net, std::vector<LPoint<>>& sampleM, std::vector<LPoint<>>& sampleB, F scan) {
        double m_Total = computeLabelTotal(sampleM.begin(), sampleM.end(), getMeasured);
        double b_Total = computeLabelTotal(sampleB.begin(), sampleB.end(), getBaseline);

        Disk maxDisk(0, 0, 0, 0);
        for (auto i = net.begin(); i != net.end() - 2; i++) {
            for (auto j = i + 1; j != net.end() - 1; j++) {
                for (auto k = j + 1; k != net.end(); k++) {
                    if (!colinear(*i, *j, *k)) {
                        double a, b, r;
                        solveCircle3(*i, *j, *k, a, b, r);
                        auto filterF = [&] (Point<> const& pt) {
                          return contains(a, b, r, pt);
                        };
                        double m_curr = computeLabelTotalF(sampleM.begin(), sampleM.end(), getMeasured,
                                        filterF);
                        double b_curr = computeLabelTotalF(sampleB.begin(), sampleB.end(), getBaseline,
                                        filterF);
                        if (scan(m_curr / m_Total, b_curr / b_Total) >= maxDisk.fValue()) {
                            maxDisk = Disk(scan(m_curr / m_Total, b_curr / b_Total), a, b, r);
                        }

                    }
                }
            }
        }
        return maxDisk;
    }

    template <typename F>
    Disk diskScan(std::vector<Point<>>& net, std::vector<Point<>>& sampleM, std::vector<Point<>>& sampleB, F scan) {
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
      std::vector <Point<>> netSampleSorted(nB, nE);

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
          auto isNotCol = [&i, &j](Point<> const& pt) {
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
          auto partitionF = [orthoX, orthoY, cX, cY](Point<> const &pt) {
              return (get<0>(pt) - cX) * orthoX + (get<1>(pt) - cY) * orthoY <= 0;
          };

          auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Point<> const &pt) {
              // If the point lines up with either of the reference
              // point then we take this to be a disk defined by only
              // the reference points.
              // We are projecting a vector created between
              //the disk center and center point between the two points.
              double a, b;
              solveCircle3(*i, *j, pt, a, b);
              return orthoX * (a - cX) + orthoY * (b - cY);
          };
          auto compF = [&orderF](Point<> const &pt1, Point<> const &pt2) {
              return orderF(pt1) < orderF(pt2);
          };

          auto asIterEnd = std::partition(msBegin, msEnd, isNotCol);
          auto bsIterEnd = std::partition(bsBegin, bsEnd, isNotCol);
          auto nIterEnd = std::partition(sortedB, sortedE, isNotCol);

          double mCount = 0;
          double bCount = 0;
          std::for_each(asIterEnd, msEnd, [i, j, &mCount](Point<> const &pt) {
              if (onLineSegment(*i, *j, pt)) {
                mCount += getMeasured(pt);
              }
          });
          std::for_each(bsIterEnd, bsEnd, [i, j, &bCount](Point<> const &pt) {
              if (onLineSegment(*i, *j, pt)) {
                bCount += getBaseline(pt);
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
          partial_counts(msBegin, mHigherIt, orderV, mCountsR, orderF, getMeasured);
          mCount = computeTotal(msBegin, mHigherIt, getMeasured) + mCount;
          partial_counts(bsBegin, bHigherIt, orderV, bCountsR, orderF, getBaseline);
          bCount = computeTotal(bsBegin, bHigherIt, getBaseline) + bCount;
          partial_counts(mHigherIt, asIterEnd, orderV, mCountsA, orderF,
                         getMeasured);
          partial_counts(bHigherIt, bsIterEnd, orderV, bCountsA, orderF,
                         getBaseline);
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
    /*
    Disk diskScanStatLabels(std::vector<LPoint<int, 2>>& net, std::vector<LPoint<int, 2>>& sampleM, std::vector<LPoint<int, 2>>& sampleB, double rho) {
        return diskScanLabels(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }
    */
    Disk diskScanStatLabels(std::vector<LPoint<>>& net, std::vector<LPoint<>>& sampleM, std::vector<LPoint<>>& sampleB, double rho) {
        return diskScanLabels(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }

    Disk diskScanSlowStatLabels(std::vector<LPoint<>>& net, std::vector<LPoint<>>& sampleM, std::vector<LPoint<>>& sampleB, double rho) {
        return diskScanSlowLabels(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }

    Disk diskScanStat(std::vector<Point<>>& net, std::vector<Point<>>& sampleM, std::vector<Point<>>& sampleB, double rho) {
        return diskScan(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }

    Disk diskScanSlowStat(std::vector<Point<>>& net, std::vector<Point<>>& sampleM, std::vector<Point<>>& sampleB, double rho) {
        return diskScanSlow(net, sampleM, sampleB, [&rho](double mr, double br){
            return kulldorff(mr, br, rho);
        });
    }
}
