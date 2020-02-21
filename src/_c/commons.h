
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
 * Reading a matrix from a file.
 * @param filename (const char *): name of the file (local path)
 * @param matrix (double **): container matrix
 * @param size (int): matrix order
 */
void read_matrix(const char * filename, double ** matrix, int size);


/*
 * Writing a matrix to a file.
 * @param filename (const char *): name of the file (local path)
 * @param matrix (double **): container matrix
 * @param size (int): matrix order
 */
void write_matrix(const char * filename, double ** matrix, int size);


/*
 * Write a vector to a file.
 * @param filename (const char *): name of the file (local path)
 * @param matrix (int *): container vector
 * @param size (int): matrix order
 */
void write_vector(const char * filename,  int * vector, int size);

#endif /* __COM */