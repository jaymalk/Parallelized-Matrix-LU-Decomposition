
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
 * @param _mat (double **): matrix
 * @param _p(int *): permutation vector
 * @param sze(int): Size of square matrices and vector
 * @param _fd(int): File descriptor
 */
void __print_permute(double ** _mat, int * _p, int sze, int _fd);


/*
 * Standard Matrix Multiplication.
 * @param _1, _2(double **): matrices involved
 * @param result(double **): result stored in this matrix
 * @param size(int): Size of the square matrices involved
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
void __matmul(double ** _1, double ** _2, double ** result, int size);

/*
 * Standard Matrix Multiplication.
 * @param original, result(double **): matrices involved
 * @param result(double **): result stored in this matrix
 * @param size(int): Size of the square matrices involved
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
double checker(double **original, double ** result, int *p, int size, int _fd);


/*
 * Reading the matrix
 */
void read_matrix(const char * filename, double ** matrix, int size);

/*
 * Write Matrix
 */
void write_matrix(const char * filename, double ** matrix, int size);

#endif /* __COM */