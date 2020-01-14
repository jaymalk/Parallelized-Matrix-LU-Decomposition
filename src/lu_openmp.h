


#ifndef __LU_OMP
#define __LU_OMP


/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (double **): input array (COPY)
 * @param l_, u_ (double **): output lu decomposition
 * @param p_ (double *): output permutation matrix
 * @param size (int): size of the arrays
 * 
 * NOTE: Copy of matrix (a_) must be passed as it gets overwritten.
 */
void __lu_decomposition(double **, double **, double **, double *, int size);


/*
 * 2D Matrix (Square), parallel initialisations.
 * @param (double ***): reference to 2D array
 * @param (int): order of the matrix
 */
void __init_2d(double ***, int);


/*
 * Complete initialisations of matrices involved (parallely).
 * @param (double ***): matrix references
 * @param (int): order of matrices
 */
void init(double ***, double ***, double ***, double **, int);


/*
 * Parallel write for matrix using openmp
 * @param (double **): matrix to write
 */
void write_matrix(double **);

#endif /* __LU_OMP */