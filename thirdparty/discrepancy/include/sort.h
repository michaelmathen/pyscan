/********************************
 * sort.h
 * Header file for sorting algorithms
 *
 * Jeff Phillips
 * jeffp@rice.edu
 * 6.4.02
 *********************************/



#ifndef SORT__H
#define SORT__H

#include <stdlib.h>
#include <stdio.h>
#include "header.h"


/***************************
 * For integer data arrays
 ***************************/

/**
 * Classic quick sort, with an integer array.
 *
 * @param data    Integer array being sorted
 * @param left    Lower bound of the part of the array being sorted
 * @param right   Upper bound of the part of the array being sorted
 */
void iquicksort (int* data, int left, int right);


/**
 * Binary search on a sorted array of integers.
 *
 * @param data    Integer array to be searched.
 * @param left    Lower bound of the part of the array being searched
 * @param right   Upper bound of the part of the array being searched
 * @param d       Data being searched for
 * @return        Index of data in array, or -1 if nonexistent.
 */
int ibinarysearch (int* data, int left, int right, int d);





/**************************************************************
 * For double data arrays paired with integer reference arrays
 **************************************************************/

/**
 * Classic quick sort, sorts the idexes of double array.
 *
 * @param data    Double array of data to be sorted over
 * @param ndx     Integer array of indexes which are actually permuted.
 * @param left    Lower bound of the part of the array being sorted
 * @param right   Upper bound of the part of the array being sorted
 */
void idquicksort (double* data, int* ndx, int left, int right);


/**
 * Binary search on a sorted array of double-integers.
 *
 * @param data    Double array to be searched.
 * @param ndx     Interger array of indexes, the part thats actually sorted
 * @param left    Lower bound of the part of the array being searched
 * @param right   Upper bound of the part of the array being searched
 * @param d       Data being searched for
 * @return        Index of data in array, or -1 if nonexistent.
 */
int idbinarysearch (double* data, int* ndx, int left, int right, double d);



#endif
