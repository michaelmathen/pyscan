/*******************************
 * interval.c
 * Source for Intervals
 * 
 * Jeff Phillips
 * jeffp@cs.duke.edu
 * 8.02.05
 ******************************/

#include "interval.h"



/**
 * Adds a point at the lowest level of an interval.
 * 
 * @param i           Interval to have pointed added.
 * @param pt          value of point
 * @param red         amount of red points
 * @param blue        amount of blue points
 * @param rw          weight of a red point
 * @param bw          weight of a blue point
 */
void add_point_ival (ival *i, double pt, double red, double blue, double rw, double bw) {

  i->d = red*rw + blue*bw;
  i->nr = (int)red;
  i->nb = (int)blue;
  i->l.d[0] = i->d;
  i->l.b[0] = pt;
  i->l.nr[0] = red;
  i->l.nb[0] = blue;
  i->r.d[0] = i->d;
  i->r.b[0] = pt;
  i->r.nr[0] = red;
  i->r.nb[0] = blue;

  if (BASE <= 0 && i->d>0) {
      i->m.d = i->d;
      i->m.l = pt;
      i->m.r = pt;
  }
}


/**
 * Adds a point to the lowest level of an interval.
 * 
 * @param i           Interval to have pointed added.
 * @param pt          value of point
 * @param red         amount of red points
 * @param blue        amount of blue points
 * @param rw          weight of a red point
 * @param bw          weight of a blue point
 */
void add_point_to_ival (ival *i, double pt, double red, double blue, double rw, double bw) {
  i->nr += (int)red;
  i->nb += (int)blue;
  i->d = i->nr*rw + i->nb*bw;
  i->l.d[0] = i->d;
  i->l.b[0] = pt;
  i->l.nr[0] = i->nr;
  i->l.nb[0] = i->nb;
  i->r.d[0] = i->d;
  i->r.b[0] = pt;
  i->r.nr[0] = i->nr;
  i->r.nb[0] = i->nb;

  if (BASE <= 0 && i->d>0) {
      i->m.d = i->d;
      i->m.l = pt;
      i->m.r = pt;
  }
}


/**
 * Merges two adjacent intervals.
 * 
 * @param n       new merged interval
 * @param l       left existing interval
 * @param r       right existing interval
 * @param rpts    to make sure we don't have > rpts-BASE
 * @param bpts    to make sure we don't have > bpts-BASE
 */
void merge_ival (ival *n, ival* l, ival* r, double rpts, double bpts) {
  int j,k;

  // overall interval
  n->d = l->d + r->d;
  n->nr = l->nr + r->nr;
  n->nb = l->nb + r->nb;
  
  // middle interval
  if (l->m.d > 0)       n->m = l->m;
  if (r->m.d > n->m.d)  n->m = r->m;  


  for (j=0; j<BSIZE; ++j) {
    for (k=0; k<BSIZE; ++k) {
      if (l->r.nr[j] + r->l.nr[k] >= BASE && 
	  l->r.nb[j] + r->l.nb[k] >= BASE) {
	if (l->r.nr[j] + r->l.nr[k] <= rpts-BASE && 
	    l->r.nb[j] + r->l.nb[k] <= bpts-BASE) {
	  if(l->r.d[j] + r->l.d[k]  >  n->m.d) {
	    n->m.d = l->r.d[j] + r->l.d[k];
	    n->m.l = l->r.b[j];
	    n->m.r = r->l.b[k];
	  }
	  else  break;
	}
	else break;
      }
    }
  }
  
  // left interval(s)
  j=0; k=0;
  do {
    if(l->l.d[j] > l->d + r->l.d[0]) {
      n->l.d[j] = l->l.d[j];
      n->l.b[j] = l->l.b[j];
      n->l.nr[j] = l->l.nr[j];
      n->l.nb[j] = l->l.nb[j];
    }
    else break;
    ++j;
  } while ((n->l.nb[j-1] < BASE && n->l.nr[j-1] < BASE) &&
	   (l->l.nb[j] != 0 || l->l.nr[j] != 0));
  
  if (j<BSIZE) {
    do {
      n->l.d[j] = l->d + r->l.d[k];
      n->l.b[j] = r->l.b[k];
      n->l.nr[j] = l->nr + r->l.nr[k];
      n->l.nb[j] = l->nb + r->l.nb[k];
      ++j; ++k;
    } while (j<BSIZE);
  }
  
  // right interval(s)
  j=0; k=0;
  do {
    if(r->r.d[j] > r->d + l->r.d[0]) {
      n->r.d[j] = r->r.d[j];
      n->r.b[j] = r->r.b[j];
      n->r.nr[j] = r->r.nr[j];
      n->r.nb[j] = r->r.nb[j];
    }
    else break;
    ++j;
  } while ((n->l.nb[j-1] < BASE && n->l.nr[j-1] < BASE) &&
	   (r->r.nb[j] != 0 || r->r.nr[j] != 0));
  
  if (j<BSIZE) {
    do {
      n->r.d[j] = r->d + l->r.d[k];
      n->r.b[j] = l->r.b[k];
      n->r.nr[j] = r->nr + l->r.nr[k];
      n->r.nb[j] = r->nb + l->r.nb[k];
      ++j; ++k;
    } while (j<BSIZE);
  }  
}


/**
 * Resets interval to have 0 ponits.
 * 
 * @param i    Interval
 */
void clear_ival (ival *i) {
  int j;
  i->d = 0;
  i->nr = 0;
  i->nb = 0;
  i->m.d = 0;
  i->m.l = i->lB;
  i->m.r = i->lB;
  for (j=0; j<BSIZE; ++j) {
    i->l.d[j] = 0;
    i->l.b[j] = i->lB;
    i->l.nr[j] = 0;
    i->l.nb[j] = 0;
    i->r.d[j] = 0;
    i->r.b[j] = i->rB;
    i->r.nr[j] = 0;
    i->r.nb[j] = 0;
  }
}


