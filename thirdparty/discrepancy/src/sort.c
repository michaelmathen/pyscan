/********************************
 * sort.c
 * Source file for sorting algorithms
 *
 * Jeff M Phillips
 * jeffp@rice.edu
 * 6.4.02
 *********************************/



#include "sort.h"


/**
 * Swaps elements of an integer array for integer quicksort.
 *
 * @param data      Interger array.
 * @param i         First index into data.
 * @param j         Second index into data.
 */
inline void iqs_swap(int* data, int i, int j) {
    int temp = data[i];
    data[i] = data[j];
    data[j] = temp;
}


/**
 * Finds a partition for integer quick sort.
 *
 * @param data         Integer array being sorted
 * @param left         Lower bound of the part of the array ebing sorted
 * @param right        Upper bound of the part of the array being sorted
 * @return             Some middle value to split over
 */
inline int iqs_partition (int* data, int left, int right) {
    int x = data[left];
    int i = left-1;
    int j = right+1;
    while (TRUE) {
        do { j--; } while (data[j] > x);
        do { i++; } while (data[i] < x);
        if (i<j) iqs_swap (data, i, j);
        else return j;
    }
}


/**
 * Classic quick sort, with an integer array.
 *
 * @param data    Integer array being sorted
 * @param left    Lower bound of the part of the array being sorted
 * @param right   Upper bound of the part of the array being sorted
 */
void iquicksort (int* data, int left, int right) {
    if (left < right) {
        int mid = iqs_partition(data, left, right);
        iquicksort (data, left, mid);
        iquicksort (data, mid+1, right);
    }
}


/**
 * Binary search on a sorted array of integers.
 *
 * @param data    Integer array to be searched.
 * @param left    Lower bound of the part of the array being searched
 * @param right   Upper bound of the part of the array being searched
 * @param d       Data being searched for
 * @return        Index of data in array, or -1 if nonexistent.
 */
int ibinarysearch (int* data, int left, int right, int d) {
    if (left < right-1) {
        int mid = (left+right)/2;
        if (d < data[mid])       return ibinarysearch (data, left, mid, d);
        else if (d > data[mid])  return ibinarysearch (data, mid+1, right, d);
        else                     return mid;
    }
    if (d == data[left])  return left;
    if (d == data[right]) return right;
    return -1;
}




/****************************************************
 * For double-integer arrays
 ****************************************************/
 

/**
 * Swaps elements of an integer array for double-integer quicksort.
 *
 * @param ndx       Integer array.
 * @param i         First index into data.
 * @param j         Second index into data.
 */
inline void idqs_swap(int* ndx, int i, int j) {
    int temp = ndx[i];
    ndx[i] = ndx[j];
    ndx[j] = temp;
}


/**
 * Finds a partition for double-integer quick sort.
 *
 * @param data         Double array.
 * @param ndx          Integer array.
 * @param left         Lower bound of the part of the array ebing sorted
 * @param right        Upper bound of the part of the array being sorted
 * @return             Some middle value to split over
 */
inline int idqs_partition (double* data, int* ndx, int left, int right) {
    double x = data[ndx[left]];
    int i = left-1;
    int j = right+1;
    while (TRUE) {
        do { j--; } while (data[ndx[j]] > x);
        do { i++; } while (data[ndx[i]] < x);
        if (i<j) idqs_swap ( ndx, i, j);
        else return j;
    }
}


/**
 * Classic quick sort, with an double-integer array.
 *
 * @param data    Double array.
 * @param ndx     Integer array.
 * @param left    Lower bound of the part of the array being sorted
 * @param right   Upper bound of the part of the array being sorted
 */
void idquicksort (double* data, int* ndx, int left, int right) {
    if (left < right) {
        int mid = idqs_partition(data, ndx, left, right);
        idquicksort (data, ndx, left, mid);
        idquicksort (data, ndx, mid+1, right);
    }
}


/**
 * Binary search on a sorted array of double-integers.
 *
 * @param data     Double array.
 * @param ndx      Integer array.
 * @param left    Lower bound of the part of the array being searched
 * @param right   Upper bound of the part of the array being searched
 * @param d       Data being searched for
 * @return        Index of data in array, or -1 if nonexistent.
 */
int idbinarysearch (double* data, int* ndx, int left, int right, double d) {
    if (left < right-1) {
        int mid = (left+right)/2;
        if (d < data[ndx[mid]])
            return idbinarysearch (data, ndx, left, mid, d);
        else if (d > data[ndx[mid]])
            return idbinarysearch (data, ndx, mid+1, right, d);
        else
            return mid;
    }
    if (d == data[ndx[left]])  return left;
    if (d == data[ndx[right]]) return right;
    return -1;
}



