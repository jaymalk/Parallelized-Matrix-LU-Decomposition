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
void __lu_decomposition(double ** a_, double **l_, double **u_, int *p_, int size) {

    // local vars
    double max; int kf;

    // local row pointers
    double *_r1, *_r2, _v, *_temp = (double *)malloc(sizeof(double)*size);

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
        swap_i(p_+k, p_+kf);
        swap_d_r(a_+k, a_+kf);
        // setting row pointers
        _r1 = l_[k], _r2 = l_[kf];
#       pragma omp parallel for
        for(int i=0; i<k; i++)
            swap_d(_r1 + i, _r2 + i);
        
        // setting values
        u_[k][k] = a_[k][k];
        double _v = u_[k][k];
        // setting row pointers
        _r1 = u_[k], _r2 = a_[k];
#       pragma omp parallel for
        for(int i=k+1; i<size; i++) {
            _temp[i] = l_[i][k] = a_[i][k]/_v;
            _r1[i] = _r2[i];
        }

        _r1 = u_[k];
#       pragma omp parallel for
        for(int i=k+1; i<size; i++) {
            _v = _temp[i];
#       pragma omp parallel for
            for(int j=k+1; j<size; j++)
                a_[i][j] -= (_v*_r1[j]);
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
void init(double *** m_, double ***mcopy, double *** l_, double *** u_, int **p_, int N) {
    // 2D init
#   pragma omp parallel sections
    {
#       pragma omp section
        __init_2d(m_, N);
#       pragma omp section
        __init_2d(mcopy, N);        
#       pragma omp section
        __init_2d(l_, N);
#       pragma omp section
        __init_2d(u_, N);
    }

    // 1D init
    (*p_) = (int *)malloc(sizeof(int)*N);

    // filling

#   pragma omp parallel for
    for(int i=0; i<N; i++) {
        // perm.matrix
    (*p_)[i] = i;
#   pragma omp parallel for
        for(int j=0; j<N; j++) {
            // matrix
            (*m_)[i][j] = drand48();
            (*mcopy)[i][j] = (*m_)[i][j];
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
    int num_threads = atoi(argv[2]);
    omp_set_num_threads(num_threads);
    double **m, **mcopy, **l, **u;
    int *p;
    double t = omp_get_wtime();
    init(&m, &mcopy, &l, &u, &p, N);
    printf("Initialization %lf\n", omp_get_wtime() - t);
    t = omp_get_wtime();
    __lu_decomposition(m, l, u, p, N);
    printf("%lf\n", omp_get_wtime() - t);
    return 0;
}
