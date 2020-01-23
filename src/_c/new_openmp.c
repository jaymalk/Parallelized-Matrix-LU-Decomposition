

#include "include.h"
#include "commons.h"

int no_of_threads = 1;

// Matrix Containers
double **a_, **l_, **u_, *p_;
const int size = 4000;


/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (double **): input array (COPY)
 * @param l_, u_ (double **): output lu decomposition
 * @param p_ (double *): output permutation matrix
 * @param size (int): size of the arrays
 * 
 * NOTE: Copy of matrix (a_) must be passed as it gets overwritten.
 */
void __lu_decomposition() {

    // local vars
    double max; int kf;

    // row pointers
    double *u_r, *a_r, *_r1, *_r2, _v;
    double *__temp = (double *)malloc(sizeof(double)*size);

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
    
        _r1 = l_[k], _r2 = l_[kf];
#       pragma omp parallel for
        for(int i=0; i<k; i++)
            swap_d(_r1 + i, _r2 + i);
        
        // setting values
        a_r = a_[k];
        u_r = u_[k];
        _v = u_r[k] = a_r[k];
#       pragma omp parallel for
        for(int i=k+1; i<size; i++) {
            __temp[i] = l_[i][k] = a_[i][k]/_v;
            u_r[i] = a_r[i];
        }
#       pragma omp parallel for
        for(int i=k+1; i<size; i++) {
#           pragma omp parallel for
            for(int j=k+1; j<size; j++)
                a_[i][j] -= (__temp[i]*u_r[j]);
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
void init() {
    // 2D init
#   pragma omp parallel sections
    {
#       pragma omp section
        __init_2d(&a_, size);
#       pragma omp section
        __init_2d(&l_, size);
#       pragma omp section
        __init_2d(&u_, size);
    }

    // 1D init
    p_ = (double *)malloc(sizeof(double)*size);

    // filling
#   pragma omp parallel
    for(int i=0; i<size; i++) {
        // perm.matrix
        p_[i] = i;
#       pragma omp parallel for
        for(int j=0; j<size; j++) {
            // matrix
            a_[i][j] = drand48();
            // u & l (conditional)
            l_[i][j] = (i>=j);
            u_[i][j] = (j>=i);
        }
    }
}



int main(int argc, char const *argv[])
{
    // num_threads
    no_of_threads = atoi(argv[1]);

    // setting threads
    omp_set_num_threads(no_of_threads);

    // init
    init();

    // timer start
    double _t = omp_get_wtime();

    // decomposition
    __lu_decomposition();


    // time set
    _t = omp_get_wtime() - _t;

    // print time
    printf("Time: %lf\n", _t);

    return 0;
}
