
#include "include.h"

#ifndef __COM
#define __COM

/*
 * Swapping double pointed values.
 * @param _1, _2 (double *)
 */
void swap_d(double *, double *);

/*
 * Swapping double pointed rows.
 * @param _1, _2 (double **)
 */
void swap_d_r(double **, double **);


/*
 * Ordered printing of a square-matrix.
 * @param _mat (double **): matrix
 * @param _sze (int): size of the matrix
 * @param _fd (file directory no.)
 * 
 * For test printing only, original write function different
 */
void _print_sq(double **, int, int);

#endif /* __COM */