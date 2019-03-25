/**
 * mtint.h  -- header for mtint.c
 *
 * Jeff M Phillips
 * jeffp@rice.edu
 * 5.21.02
 */

#include <stdio.h>
#include <math.h>
#include "header.h"

/**
 * Initializes the RNG with a seed.
 *
 * @param seed The seed.
 */
void sgenrand (unsigned long seed);


/**
 * Auto seeds the RNG (maybe with random stuff from hardware.
 */
void sgenrand_auto(void);

/**
 * Seeds RNG and then stores seed.
 *
 * @param seed  The seed.
 * @param sel   If true, pump through RNG once.
 * @param fp    File to save seeds in.
 */
void sgenrand_sel(unsigned long seed, int sel, FILE *fp);


/**
 * Seeds RNF with an array of unsigned long integers.
 *
 * @param seed_array The array of unsigned long integers.
 */
void lsgenrand(unsigned long seed_array[]);


/**
 * Generates a random unsigned long.
 *
 * @return random number.
 */
unsigned long genrand (void);

/**
 * Dumps the bits of an input.
 *
 * @param x input.
 */
void dumpbits(unsigned long x);


/**
 * Generates a random number out of a uniform distribution.
 *
 * @return random number
 */
double genuniform (void);


/**
 * Generates a random number out of a N(0,1) normal distribution
 * 
 * @return a normal random number
 */
double gennormal (void);
