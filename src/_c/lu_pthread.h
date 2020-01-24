

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
 * @param (int): order of matrices
 */
void init(int);

//Pointer function for pthread to swap values of l[k] and l[k']
void *swap_l(void *);

//Pointer function for pthread to compute values of l and u
void *lu(void *);

//Pointer function for pthread to compute values of matrix a
void *mlu(void *);

//Pointer function for pthread to randomly initialize the matrices
void *threaded_init(void *);

//Pointer function for pthread to allocate memory to each row
void *threaded_alloc(void *);

#endif /* __LU_PTHREADS */