

#include "include.h"
#include "commons.h"

#include "lu_openmp.h"



/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (double **): input array (COPY)
 * @param l_, u_ (double **): output lu decomposition
 * @param p_ (double *): output permutation matrix
 * @param size (int): size of the arrays
 * 
 * NOTE: Copy of matrix (a_) must be passed as it gets overwritten.
 */
void __lu_decomposition(double ** a_, double **l_, double **u_, double *p_, int size) {

    // local vars
    double max; int kf;

    // basic loop (row-wise)
    for(int k = 0; k < size; k++) {
        max = 0;
        for(int i=k; i<size; i++) {
            if(max < fabs(a_[i][k])) {
                max = fabs(a_[i][k]);
                kf = i;
            }
        }
        if(max == 0) {
            fprintf(__stderrp, "Singular Matrix...\nExiting\n");
            exit(2);
        }

        // swapping
        swap_d(p_+k, p_+kf);
        swap_d_r(a_+k, a_+kf);
    
#       pragma omp parallel for
        for(int i=0; i<k; i++)
            swap_d(l_[k]+i, l_[kf]+i);
        
        // setting values
        u_[k][k] = a_[k][k];
#       pragma omp parallel for
        for(int i=k+1; i<size; i++) {
            l_[i][k] = a_[i][k]/u_[k][k];
            u_[k][i] = a_[k][i];
        }
#       pragma omp parallel for
        for(int i=k+1; i<size; i++) {
#           pragma omp parallel for
            for(int j=k+1; j<size; j++)
                a_[i][j] -= (l_[i][k]*u_[k][j]);
        }
    }
}


/*
 * 2D Matrix (Square), parallel initialisations.
 * @param (double ***): reference to 2D array
 * @param (int): order of the matrix
 */
void __init_2d(double *** _m, int _sze) {
    (*_m) = (double **)malloc(sizeof(double *)*_sze);
#   pragma omp parallel for
    for(int i=0; i<_sze; i++)
        (*_m)[i] = (double *)malloc(sizeof(double)*_sze);
}


/*
 * Complete initialisations of matrices involved (parallely).
 * @param (double ***): matrix references
 * @param (int): order of matrices
 */
void init(double *** m_, double *** l_, double *** u_, double **p_, int N) {
    // 2D init
#   pragma omp parallel sections
    {
#       pragma omp section
        __init_2d(m_, N);
#       pragma omp section
        __init_2d(l_, N);
#       pragma omp section
        __init_2d(u_, N);
    }

    // 1D init
    (*p_) = (double *)malloc(sizeof(double)*N);

    // filling
#   pragma omp parallel for
    for(int i=0; i<N; i++) {
        // perm.matrix
        (*p_)[i] = i;
#       pragma omp parallel for
        for(int j=0; j<N; j++) {
            // matrix
            (*m_)[i][j] = drand48();
            // u & l (conditional)
            if(i == j) {
                (*l_)[i][j] = 1;
                (*u_)[i][j] = 1;
            }
            else if (i > j) {
                (*l_)[i][j] = 1;
            }
            else {
                (*u_)[i][j] = 1;
            }
        }
    }
}



int main(int argc, char const *argv[])
{
    int N = atoi(argv[1]);
    double **m, **l, **u, *p;
    init(&m, &l, &u, &p, N);
    __lu_decomposition(m, l, u, p, N);
    _print_sq(l, N, 2);
    return 0;
}
