//
// Created by mmath on 9/25/17.
//
#include <algorithm>
#include <unordered_map>
#include <ostream>

#include <set>
#include <deque>
#include <unordered_set>
#include <assert.h>

#include "SparseGrid.hpp"
#include "Statistics.hpp"
#include "Point.hpp"
#include "DiskScan.hpp"

namespace pyscan {

    bool colinear(Pt2 const& pt1,
                  Pt2 const& pt2,
                  Pt2 const& pt3){
        double x1, x2, x3, y1, y2, y3;
        getLoc(pt1, x1, y1);
        getLoc(pt2, x2, y2);
        getLoc(pt3, x3, y3);

        double
                a11 = x2 - x1, a12 = y2 - y1,
                a21 = x2 - x3, a22 = y2 - y3;
        return (a11 * a22 - a12 * a21 == 0);
    }


    bool onLineSegment(Pt2 const& pt1,
                       Pt2 const& pt2,
                       Pt2 const& pt3) {
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

    void solveCircle3(Pt2 const& pt1, Pt2 const& pt2, Pt2 const& pt3,
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

    void solveCircle3(Pt2 const& pt1, Pt2 const& pt2, Pt2 const& pt3,
                     double &a, double &b, double &r) {
        double x1, x2, x3, y1, y2, y3;
        getLoc(pt1, x1, y1);
        getLoc(pt2, x2, y2);
        getLoc(pt3, x3, y3);
        // Setup a matrix equation of the form Ax = b

        // A
        double a11 = x2 - x1,
                a12 = y2 - y1,
                a21 = x2 - x3,
                a22 = y2 - y3;
        // b
        double b1 = (y2 * y2 + x2 * x2 - y1 * y1 - x1 * x1) / 2,
                b2 = (y2 * y2 + x2 * x2 - y3 * y3 - x3 * x3) / 2;

        double detA = a11 * a22 - a12 * a21;
        // inverse of A
        double ai11 = a22 / detA,
                ai12 = -a12 / detA,
                ai21 = -a21 / detA,
                ai22 = a11 / detA;

        // A^{-1} b = x
        a = ai11 * b1 + ai12 * b2;
        b = ai21 * b1 + ai22 * b2;
        r = sqrt((x1 - a) * (x1 - a) + (y1 - b) * (y1 - b));
    }

    void findPerpVect(Pt2 const& p1, Pt2 const& p2, double* u, double* v) {
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





    template<typename T, typename F>
    inline void partial_counts(T begin, T end,
                               std::vector<double> const& partitions,
                               std::vector<double>& counts, F orderF) {
        //Partitions based on the break points.
        for (; begin != end; begin++) {
            auto lb = std::lower_bound(partitions.begin(), partitions.end(),
                                  orderF(*begin));
            counts[lb - partitions.begin()] += getWeight(*begin);
        }
    }


    template<typename T, typename F>
    inline void partial_counts_label(T begin, T end,
                               std::vector<double> const& partitions,
                               std::vector<crescent_t >& counts, F orderF) {
        //Partitions based on the break points.
        for (; begin != end; begin++) {
            auto lb = std::lower_bound(partitions.begin(), partitions.end(),
                                  orderF(*begin));
            counts[lb - partitions.begin()].emplace_back(begin->label, getWeight(*begin));
        }
    }

    double updateCounts(std::unordered_map<size_t, size_t>& curr_counts,
                        crescent_t const& adding, crescent_t const& removing) {
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

    std::ostream& operator<<(std::ostream& os, std::unordered_map<size_t, size_t> const& items) {
      for (auto ix : items) {
        os << ix.first << ":" << ix.second << ", ";
      }
      return os;
    }



    std::ostream& operator<<(std::ostream& os, crescent_t const& items) {
      for (auto ix : items) {
        os << ix.label  << ", ";
      }
      return os;
    }


    std::tuple<Disk, double> diskScanLabels(point_list& net,
                        lpoint_list& sampleM,
                        lpoint_list& sampleB,
                        std::function<double(double, double)> const& scan) {
        /*
         * Each point is defined with a label. For each division we store the set of
         * labels.
         */
        //Calculate the total measured and baseline value.
        auto msBegin = sampleM.begin();
        auto bsBegin = sampleB.begin();
        auto msEnd = sampleM.end();
        auto bsEnd = sampleB.end();
        double m_Total = computeLabelTotal(msBegin, msEnd);
        double b_Total = computeLabelTotal(bsBegin, bsEnd);

        auto nB = net.begin();
        auto nE = net.end();
        std::vector <Point<>> netSampleSorted(nB, nE);

        Disk currMax;
        double maxStat = 0;

        auto sortedB = netSampleSorted.begin();
        auto sortedEOuter = netSampleSorted.end();
        for (auto i = nB; i != nE - 2; i++) {
            sortedEOuter = std::remove(sortedB, sortedEOuter, *i);
            auto sortedE = sortedEOuter;
            for (auto j = i + 1; j != nE - 1; j++) {
                auto el = std::find(sortedB, sortedE, *j);
                if (el == sortedE) {
                    continue;
                }
                std::swap(*el, *(sortedE - 1));
                sortedE = sortedE - 1;
                //Create a vector between the two points
                double orthoX, orthoY;
                findPerpVect(*i, *j, &orthoX, &orthoY);
                double cX = (getX(*i) + getX(*j)) / 2.0;
                double cY = (getY(*i) + getY(*j)) / 2.0;
                auto isNotCol = [&i, &j](Point<> const& pt) {
                    return !colinear(*i, *j, pt);
                };
                // Partition these into a set of adding points and removing points
                auto partitionF = [orthoX, orthoY, cX, cY](Point<> const &pt) {
                    return (getX(pt) - cX) * orthoX + (getY(pt) - cY) * orthoY <= 0;
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

                // will have two sets an added set and a current set.
                // added set is for stuff that will never leave. The current set on the other hand
                // will correspond to every trajectory that we overlap.
                auto onSegment = [&](Point<> const& pt) {
                    return onLineSegment(*i, *j, pt);
                };

                //Counts the points that lie on the line segment between i and j.
                //These points are colinear so have been removed from the scan.
                std::unordered_map<size_t, size_t> m_curr_set;
                std::unordered_map<size_t, size_t> b_curr_set;
                double m_count = computeLabelTotalF(asIterEnd, msEnd, m_curr_set, onSegment);
                double b_count = computeLabelTotalF(bsIterEnd, bsEnd, b_curr_set, onSegment);

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
                partial_counts_label(msBegin, mHigherIt, orderV, mCountsR, orderF);
                m_count += computeLabelTotal(msBegin, mHigherIt, m_curr_set);
                partial_counts_label(bsBegin, bHigherIt, orderV, bCountsR, orderF);
                b_count += computeLabelTotal(bsBegin, bHigherIt, b_curr_set);
                partial_counts_label(mHigherIt, asIterEnd, orderV, mCountsA, orderF);
                partial_counts_label(bHigherIt, bsIterEnd, orderV, bCountsA, orderF);
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
                        Disk currDisk(a, b, r);
                        currMax = currDisk;
                        maxStat = newStat;
                    }
                }
            }
        }
        return std::make_tuple(currMax, maxStat);
    }

    bool contains(double a, double b, double r, Pt2 const& pt) {
        double la = (getX(pt) - a);
        double lb = (getY(pt) - b);
        return la * la + lb * lb <= r * r;
    }

    auto evaluateRegion(pyscan::wpoint_list const& m_pts, pyscan::wpoint_list const& b_pts,
                        pyscan::Disk const& disk, std::function<double(double, double)> const& scan) -> double {

        double m_val = 0 ,m_total = 0, b_val = 0, b_total = 0;
        std::for_each(m_pts.begin(), m_pts.end(), [&](auto const& pt){
            if (disk.contains(pt)) {
                m_val += pyscan::getWeight(pt);
            }
            m_total += pyscan::getWeight(pt);
        });
        std::for_each(b_pts.begin(), b_pts.end(), [&](auto const& pt){
            if (disk.contains(pt)) {
                b_val += pyscan::getWeight(pt);
            }
            b_total += pyscan::getWeight(pt);
        });
        return scan(m_val / m_total, b_val / b_total);
    }

    auto evaluateRegion(pyscan::lpoint_list const& m_pts, pyscan::lpoint_list const& b_pts,
                        pyscan::Disk const& disk, std::function<double(double, double)> const& scan)-> double {

        double m_total = computeLabelTotal(m_pts.begin(), m_pts.end());
        double m_val = computeLabelTotalF(m_pts.begin(), m_pts.end(), [&](auto const& pt) {
             return disk.contains(pt);
        });

        double b_total = computeLabelTotal(b_pts.begin(), b_pts.end());
        double b_val = computeLabelTotalF(b_pts.begin(), b_pts.end(), [&](auto const& pt) {
            return disk.contains(pt);
        });
        return scan(m_val / m_total, b_val / b_total);
    }

    std::tuple<Disk, double> diskScanSlow(point_list const& net, wpoint_list const& sampleM, wpoint_list const& sampleB,
                      std::function<double(double, double)> const& scan) {
        double m_Total = computeTotal(sampleM.begin(), sampleM.end());
        double b_Total = computeTotal(sampleB.begin(), sampleB.end());

        Disk maxDisk;
        double max_stat = 0;
        for (auto i = net.begin(); i != net.end() - 2; i++) {
            for (auto j = i + 1; j != net.end() - 1; j++) {
                for (auto k = j + 1; k != net.end(); k++) {
                    if (!colinear(*i, *j, *k)) {
                        double a, b, r;
                        solveCircle3(*i, *j, *k, a, b, r);
                        Disk disk(a, b, r);
                        double m_curr = 0;
                        double b_curr = 0;
                        for (auto& p : sampleM) {
                            if (contains(a, b, r, p)) {
                                m_curr += getWeight(p);
                            }
                        }

                        for (auto& p : sampleB) {
                            if (contains(a, b, r, p)) {
                                b_curr += getWeight(p);
                            }
                        }

                        if (scan(m_curr / m_Total, b_curr / b_Total) >= max_stat) {
                            maxDisk = Disk(a, b, r);
                            max_stat = scan(m_curr / m_Total, b_curr / b_Total);
                        }
                    }
                }
            }
        }
        return std::make_tuple(maxDisk, max_stat);
    }

    std::tuple<Disk, double> diskScanSlowLabels(point_list& net, lpoint_list& sampleM, lpoint_list& sampleB,
                            std::function<double(double, double)> const& scan) {
        double m_Total = computeLabelTotal(sampleM.begin(), sampleM.end());
        double b_Total = computeLabelTotal(sampleB.begin(), sampleB.end());

        Disk maxDisk(0, 0, 0);
        double max_scan = 0;
        for (auto i = net.begin(); i != net.end() - 2; i++) {
            for (auto j = i + 1; j != net.end() - 1; j++) {
                for (auto k = j + 1; k != net.end(); k++) {
                    if (!colinear(*i, *j, *k)) {
                        double a, b, r;
                        solveCircle3(*i, *j, *k, a, b, r);
                        auto filterF = [&] (Pt2 const& pt) {
                          return contains(a, b, r, pt);
                        };
                        double m_curr = computeLabelTotalF(sampleM.begin(), sampleM.end(),
                                        filterF);
                        double b_curr = computeLabelTotalF(sampleB.begin(), sampleB.end(),
                                        filterF);
                        if (scan(m_curr / m_Total, b_curr / b_Total) >= max_scan) {
                            maxDisk = Disk( a, b, r);
                            max_scan = scan(m_curr / m_Total, b_curr / b_Total);
                        }
                    }
                }
            }
        }
        return std::make_tuple(maxDisk, max_scan);
    }

    std::tuple<Disk, double> diskScan(point_list& net,
                  wpoint_list& sampleM,
                  wpoint_list& sampleB,
                  std::function<double(double, double)> const& scan) {
      //Calculate the total measured and baseline value.
      double m_Total = 0;
      double b_Total = 0;
      auto msBegin = sampleM.begin();
      auto bsBegin = sampleB.begin();
      auto msEnd = sampleM.end();
      auto bsEnd = sampleB.end();
      std::for_each(msBegin, msEnd, [&](WPoint<> const& pt) {
          m_Total += getWeight(pt);
      });

      std::for_each(bsBegin, bsEnd, [&](WPoint<> const& pt) {
          b_Total += getWeight(pt);
      });

      auto nB = net.begin();
      auto nE = net.end();
      std::vector <Pt2> netSampleSorted(nB, nE);

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
          if (el == sortedE) {
              continue;
          }
          std::swap(*el, *(sortedE - 1));
          sortedE = sortedE - 1;
          //Create a vector between the two points
          double orthoX, orthoY;
          findPerpVect(*i, *j, &orthoX, &orthoY);
          double cX = (getX(*i) + getX(*j)) / 2.0;
          double cY = (getY(*i) + getY(*j)) / 2.0;
          auto isNotCol = [&i, &j](Pt2 const& pt) {
            return !colinear(*i, *j, pt);
          };
          // Partition these into a set of adding points and removing points
          auto partitionF = [orthoX, orthoY, cX, cY](Pt2 const &pt) {
              return (getX(pt) - cX) * orthoX + (getY(pt) - cY) * orthoY <= 0;
          };

          auto orderF = [orthoX, orthoY, &i, &j, cX, cY](Pt2 const &pt) {
              // If the point lines up with either of the reference
              // point then we take this to be a disk defined by only
              // the reference points.
              // We are projecting a vector created between
              //the disk center and center point between the two points.
              double a, b;
              solveCircle3(*i, *j, pt, a, b);
              return orthoX * (a - cX) + orthoY * (b - cY);
          };
          auto compF = [&orderF](Pt2 const &pt1, Pt2 const &pt2) {
              return orderF(pt1) < orderF(pt2);
          };

          auto asIterEnd = std::partition(msBegin, msEnd, isNotCol);
          auto bsIterEnd = std::partition(bsBegin, bsEnd, isNotCol);
          auto nIterEnd = std::partition(sortedB, sortedE, isNotCol);

          double mCount = 0;
          double bCount = 0;
          std::for_each(asIterEnd, msEnd, [i, j, &mCount](WPoint<> const &pt) {
              if (onLineSegment(*i, *j, pt)) {
                mCount += getWeight(pt);
              }
          });
          std::for_each(bsIterEnd, bsEnd, [i, j, &bCount](WPoint<> const &pt) {
              if (onLineSegment(*i, *j, pt)) {
                bCount += getWeight(pt);
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
          partial_counts(msBegin, mHigherIt, orderV, mCountsR, orderF);
          mCount = computeTotal(msBegin, mHigherIt) + mCount;
          partial_counts(bsBegin, bHigherIt, orderV, bCountsR, orderF);
          bCount = computeTotal(bsBegin, bHigherIt) + bCount;
          partial_counts(mHigherIt, asIterEnd, orderV, mCountsA, orderF);
          partial_counts(bHigherIt, bsIterEnd, orderV, bCountsA, orderF);
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
                Disk currDisk(a, b, r);
                currMax = currDisk;
                maxStat = newStat;
            }
          }
        }
      }
      return std::make_tuple(currMax, maxStat);
    }


    template<typename T, typename TW>
    std::tuple<Disk, double> disk_scan_restricted(Point<> const& p1,
                                                Point<> const& p2,
                                                T& net,
                                                TW& sampleM,
                                                TW& sampleB,
                                                double min_dist,
                                                double max_dist,
                                                double m_Total,
                                                double b_Total,
                                                std::function<double(double, double)> const& scan) {


        Disk currMax;
        double maxStat = 0;

        auto nB = net.begin(), nE = net.end();
        auto msBegin = sampleM.begin(), msEnd = sampleM.end();
        auto bsBegin = sampleB.begin(), bsEnd = sampleB.end();

        //Create a vector between the two points
        double orthoX, orthoY;
        findPerpVect(p1, p2, &orthoX, &orthoY);
        double cX = (getX(p1) + getX(p2)) / 2.0;
        double cY = (getY(p1) + getY(p2)) / 2.0;
        auto isNotCol = [&p1, &p2](Pt2 const& pt) {
            return !colinear(p1, p2, pt);
        };
        // Partition these into a set of adding points and removing points
        auto partitionF = [orthoX, orthoY, cX, cY](Pt2 const &pt) {
            return (getX(pt) - cX) * orthoX + (getY(pt) - cY) * orthoY <= 0;
        };

        auto orderF = [orthoX, orthoY, p1, p2, cX, cY](Pt2 const &pt) {
            // If the point lines up with either of the reference
            // point then we take this to be a disk defined by only
            // the reference points.
            // We are projecting a vector created between
            //the disk center and center point between the two points.
            double a, b;
            solveCircle3(p1, p2, pt, a, b);
            return orthoX * (a - cX) + orthoY * (b - cY);
        };
        auto compF = [&orderF](Pt2 const &pt1, Pt2 const &pt2) {
            return orderF(pt1) < orderF(pt2);
        };

        auto asIterEnd = std::partition(msBegin, msEnd, isNotCol);
        auto bsIterEnd = std::partition(bsBegin, bsEnd, isNotCol);
        nE = std::partition(nB, nE, isNotCol);

        auto tooBigTooSmall = [&] (Point<> const& pt) {
            double a, b, r;
            solveCircle3(p1, p2, pt, a, b, r);
            return min_dist < r && r <= max_dist;
        };

        nE = std::partition(nB, nE, tooBigTooSmall);

        double mCount = 0;
        double bCount = 0;
        std::for_each(asIterEnd, msEnd, [p1, p2, &mCount](WPoint<> const &pt) {
            if (onLineSegment(p1, p2, pt)) {
                mCount += getWeight(pt);
            }
        });
        std::for_each(bsIterEnd, bsEnd, [p1, p2, &bCount](WPoint<> const &pt) {
            if (onLineSegment(p1, p2, pt)) {
                bCount += getWeight(pt);
            }
        });

        std::sort(nB, nE, compF);
        std::vector<double> orderV;
        orderV.reserve(nE - nB);
        for (auto b = nB; b != nE; b++) {
            orderV.push_back(orderF(*b));
        }
        auto mHigherIt = std::partition(msBegin, asIterEnd, partitionF);
        auto bHigherIt = std::partition(bsBegin, bsIterEnd, partitionF);

        std::vector<double> mCountsA(nE - nB, 0);
        std::vector<double> bCountsA(nE - nB, 0);
        std::vector<double> mCountsR(nE - nB, 0);
        std::vector<double> bCountsR(nE - nB, 0);
        std::fill(mCountsR.begin(), mCountsR.end(), 0);
        std::fill(bCountsR.begin(), bCountsR.end(), 0);
        std::fill(mCountsA.begin(), mCountsA.end(), 0);
        std::fill(bCountsA.begin(), bCountsA.end(), 0);
        /*Probably most of the time is spent here*/
        partial_counts(msBegin, mHigherIt, orderV, mCountsR, orderF);
        mCount = computeTotal(msBegin, mHigherIt) + mCount;
        partial_counts(bsBegin, bHigherIt, orderV, bCountsR, orderF);
        bCount = computeTotal(bsBegin, bHigherIt) + bCount;
        partial_counts(mHigherIt, asIterEnd, orderV, mCountsA, orderF);
        partial_counts(bHigherIt, bsIterEnd, orderV, bCountsA, orderF);
        /*----------------------------------------------*/
        //Now scan over the counts.
        auto size = nE - nB;
        for (int k = 0; k < size; k++) {
            mCount += mCountsA[k] - mCountsR[k];
            bCount += bCountsA[k] - bCountsR[k];
            double m_hat = mCount / m_Total;
            double b_hat = bCount / b_Total;
            double newStat = scan(m_hat, b_hat);
            if (maxStat <= newStat) {
                double a, b, r;
                solveCircle3(p1, p2, *(nB + k), a, b, r);
                Disk currDisk(a, b, r);
                currMax = currDisk;
                maxStat = newStat;
            }
        }
        return std::make_tuple(currMax, maxStat);
    }


    template <typename P, typename T>
    void insert_queue(std::deque<P> &items, T b_it, T e_it) {
        for (auto b = b_it; b != e_it; b++) {
            items.push_front(*b);
        }
    }

    template <typename P, typename T>
    void remove_queue(std::deque<P> &items, T b_it, T e_it) {
        for (auto b = b_it; b != e_it; b++) {
            items.pop_back();
        }
    }

    std::tuple<Disk, double> disk_scan_scale(point_list& net,
                                            wpoint_list& sampleM,
                                            wpoint_list& sampleB,
                                            uint32_t grid_r,
                                            std::function<double(double, double)> const& scan) {

            //Calculate the total measured and baseline value.
            double m_Total = std::accumulate(sampleM.begin(), sampleM.end(), 0.0, [&](double w, WPoint<> const& p){
                return w + getWeight(p);
            });

            double b_Total = std::accumulate(sampleB.begin(), sampleB.end(), 0.0, [&](double w, WPoint<> const& pt) {
                return w + getWeight(pt);
            });

            SparseGrid<Point<>> grid_net(net, grid_r);
            SparseGrid<WPoint<>> grid_sample_m(sampleM, grid_r);
            SparseGrid<WPoint<>> grid_sample_b(sampleB, grid_r);

            using net_it = typename SparseGrid<Point<>>::pt_it;
            using pt_it = typename SparseGrid<WPoint<>>::pt_it;
            Disk currMax;
            double maxStat = 0;


            for (uint32_t i = 0; i < grid_r; i++) { // through y

                std::deque<Point<>> net_chunk;
                std::deque<WPoint<>> m_sample_chunk;
                std::deque<WPoint<>> b_sample_chunk;

                for (int k = i - 3; k <= i + 1; k++) {// Through y, but ignore the last part
                    for (int l = 0; l <= 2; l++) { // through x
                        //Rotate this queue
                        net_it chunk_begin, chunk_end;
                        std::tie(chunk_begin, chunk_end) = grid_net(k, l);
                        insert_queue(net_chunk, chunk_begin, chunk_end);

                        pt_it wchunk_begin, wchunk_end;
                        std::tie(wchunk_begin, wchunk_end) = grid_sample_m(k, l);
                        insert_queue(m_sample_chunk, wchunk_begin, wchunk_end);

                        std::tie(wchunk_begin, wchunk_end) = grid_sample_b(k, l);
                        insert_queue(b_sample_chunk, wchunk_begin, wchunk_end);
                    }
                }
                for (int j = 0; j < grid_r; j++) { // through x

                    // extend this vector with points with this entire block
                    for (int k = i - 2; k <= i + 2; i++) {

                        //Rotate this queue
                        net_it chunk_begin, chunk_end;
                        std::tie(chunk_begin, chunk_end) = grid_net(j + 2, k);
                        insert_queue(net_chunk, chunk_begin, chunk_end);

                        pt_it wchunk_begin, wchunk_end;
                        std::tie(wchunk_begin, wchunk_end) = grid_sample_m(j + 2, k);
                        insert_queue(m_sample_chunk, wchunk_begin, wchunk_end);

                        std::tie(wchunk_begin, wchunk_end) = grid_sample_b(j + 2, k);
                        insert_queue(b_sample_chunk, wchunk_begin, wchunk_end);


                        std::tie(chunk_begin, chunk_end) = grid_net(j - 3, k);
                        remove_queue(net_chunk, chunk_begin, chunk_end);

                        std::tie(wchunk_begin, wchunk_end) = grid_sample_m(j - 3, k);
                        remove_queue(m_sample_chunk, wchunk_begin, wchunk_end);

                        std::tie(wchunk_begin, wchunk_end) = grid_sample_b(j - 3, k);
                        remove_queue(b_sample_chunk, wchunk_begin, wchunk_end);
                    }

                    net_it pt_begin, pt_end;
                    std::tie(pt_begin, pt_end) = grid_net(i, j);
                    for (auto pt1 = pt_begin; pt1 != pt_end - 2; pt1++) {
                        for (auto pt2 = pt1 + 1; pt2 != pt_end - 1; pt2++) {
                            Disk d1;
                            double possible_max;
                            std::tie(d1, possible_max) = disk_scan_restricted(*pt1, *pt2,
                                    net_chunk,
                                    m_sample_chunk,
                                    b_sample_chunk,
                                    grid_net.get_resolution(),
                                    2 * grid_net.get_resolution(),
                                    m_Total,
                                    b_Total,
                                    scan);
                            if (possible_max > maxStat) {
                                currMax = d1;
                                maxStat = possible_max;
                            }
                        }
                    }
                }
            }
            return std::make_tuple(currMax, maxStat);
    }
}
