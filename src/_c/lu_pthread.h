

#ifndef __LU_PTHREADS
#define __LU_PTHREADS

#include <pthread.h>
#include <omp.h>
/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (double **): input array (COPY)
 * @param l_, u_ (double **): output lu decomposition
 * @param p_ (double *): output permutation matrix
 * @param size (int): size of the arrays
 * 
 * NOTE: Copy of matrix (a_) must be passed as it gets overwritten.
 */
void __lu_decomposition(int);

/*
 * 2D Matrix (Square), parallel initialisations.
 * @param (double ***): reference to 2D array
 * @param (int): order of the matrix
 */
void *__init_2d(void *);


/*
 * Complete initialisations of matrices involved (parallely).
 * @param (double ***): matrix references
 * @param (int): order of matrices
 */
void init(int);


/*
 * Parallel write for matrix using pthreads
 * @param (double **): matrix to write
 */

void *swap_l(void *);
void *lu(void *);
void *mlu(void *);
void *threaded_alloc(void *);
void *threaded_init(void *);

#endif /* __LU_PTHREADS */