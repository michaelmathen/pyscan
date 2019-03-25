/******************************
 * interval.h
 * Data structure for maintaining discrepancy intervals
 *
 * Jeff Phillips
 * jeffp@cs.duke.edu
 * 8.02.05
 ******************************/

#ifndef INTERVAL__H
#define INTERVAL__H


//#define BASE 3
#define BSIZE BASE*2

#include "header.h"

///////////// STRUCTURE //////////////////

struct maximum_interval_t {
  double l, r;        // left and right boundary of interval
  double d;           // discrepancy of interval
};
typedef struct maximum_interval_t m_ival;

struct boundary_interval_t {
  double b[BSIZE+1];     // variable boundary
  double d[BSIZE+1];     // discrepancy at each boundary
  double nr[BSIZE+1];    // number of red points in [fb,vb[i]]
  double nb[BSIZE+1];    // number of blue points in [fb,vb[i]]
};
typedef struct boundary_interval_t b_ival;

struct interval_t {
  double d;               // discrepancy of interval
  double lB,rB;           // left and right boundary
  m_ival m;               // maximum interval
  b_ival l,r;             // left, right interval
  int nr, nb;             // number of red and blue points in the interval
};
typedef struct interval_t ival;


//////////// FUCNTIONS ///////////////////

void add_point_ival (ival *i, double pt, double red, double blue, double rw, double bw);
void add_point_to_ival (ival *i, double pt, double red, double blue, double rw, double bw);
void merge_ival (ival *n, ival* l, ival* r, double rpts, double bpts);
void clear_ival (ival *i);

#endif
