
#include "include.h"

#ifndef __COM
#define __COM

/*
 * Swapping double pointed values.
 * @param _1, _2 (double *)
 */
void swap_d(double *, double *);

/*
 * Swapping int pointed values.
 * @param _1, _2 (int *)
 */
void swap_i(int * _1, int * _2); 

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

/*
 * Ordered printing of a square-matrix, keeping in mind the permutation.
 * @param _mat (vector<vector<double>>): matrix
 */
void __print_permute(double ** _mat, int * _p, int sze, int _fd);


/*
 * Standard Matrix Multiplication.
 * @param (std::vector<std::vector<double>>&, std::vector<std::vector<double>>&): matrices involved
 * 
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
void __matmul(double ** _1, double ** _2, double ** result, int size);

double checker(double **original, double ** result, int *p, int size, int);
#endif /* __COM */