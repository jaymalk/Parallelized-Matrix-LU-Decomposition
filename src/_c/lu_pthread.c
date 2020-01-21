#include "commons.h"
#include "lu_pthread.h"

int no_of_threads = 1;
double **m, **l, **u, *p;

struct init_arg_struct {
    double ** mat;
    int size;
};

/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (double **): input array (COPY)
 * @param l_, u_ (double **): output lu decomposition
 * @param p_ (double *): output permutation matrix
 * @param size (int): size of the arrays
 * 
 * NOTE: Copy of matrix (a_) must be passed as it gets overwritten.
 */
void *__lu_decomposition(int size) {

    // local vars
    double max; int kf;

    // basic loop (row-wise)
    for(int k = 0; k < size; k++) {
        max = 0;
        for(int i=k; i<size; i++) {
            if(max < fabs(m[i][k])) {
                max = fabs(m[i][k]);
                kf = i;
            }
        }
        if(max == 0) {
            fprintf(__stderrp, "Singular Matrix...\nExiting\n");
            exit(2);
        }

        // swapping
        swap_d(p+k, p+kf);
        swap_d_r(m+k, m+kf);
    
#       pragma omp parallel for num_threads(no_of_threads)
        for(int i=0; i<k; i++)
            swap_d(l[k]+i, l[kf]+i);
        
        // setting values
        u[k][k] = m[k][k];
#       pragma omp parallel for num_threads(no_of_threads)
        for(int i=k+1; i<size; i++) {
            l[i][k] = m[i][k]/u[k][k];
            u[k][i] = m[k][i];
        }
#       pragma omp parallel for num_threads(no_of_threads)
        for(int i=k+1; i<size; i++) {
#           pragma omp parallel for
            for(int j=k+1; j<size; j++)
                m[i][j] -= (l[i][k]*u[k][j]);
        }
    }
}


/*
 * 2D Matrix (Square), parallel initialisations.
 * @param (double ***): reference to 2D array
 * @param (int): order of the matrix
 */
void *__init_2d(void *arguments) {

    struct init_arg_struct *args = (struct init_arg_struct *)arguments;
    double ** _m = args->mat;
    int _sze = args->size;
    (_m) = (double **)malloc(sizeof(double *)*_sze);
#   pragma omp parallel for num_threads(no_of_threads)
    for(int i=0; i<_sze; i++)
        (_m)[i] = (double *)malloc(sizeof(double)*_sze);
}


/*
 * Complete initialisations of matrices involved (parallely).
 * @param (double ***): matrix references
 * @param (int): order of matrices
 */
void init(int N) {
    // 2D init
    pthread_t init_thread_id[3];
    {
        struct init_arg_struct args;
        args.size = N;
        args.mat = m;
        pthread_create( &init_thread_id[0], NULL, &__init_2d, (void *)&args);
        
        args.mat = l;
        pthread_create( &init_thread_id[1], NULL, &__init_2d, (void *)&args);
        
        args.mat = u;
        pthread_create( &init_thread_id[2], NULL, &__init_2d, (void *)&args);
        
    }

    // 1D init
    (p) = (double *)malloc(sizeof(double)*N);

    // filling
#   pragma omp parallel for num_threads(no_of_threads)
    for(int i=0; i<N; i++) {
        // perm.matrix
        (p)[i] = i;
#       pragma omp parallel for
        for(int j=0; j<N; j++) {
            // matrix
            (m)[i][j] = drand48();
            // u & l (conditional)
            if(i == j) {
                (l)[i][j] = 1;
                (u)[i][j] = 1;
            }
            else if (i > j) {
                (l)[i][j] = 1;
            }
            else {
                (u)[i][j] = 1;
            }
        }
    }
}



int main(int argc, char const *argv[])
{
    int N = atoi(argv[1]);
    no_of_threads = atoi(argv[2]);
    
    double t = omp_get_wtime();
    init(N);
    printf("Initialization %lf\n", omp_get_wtime() - t);
    
    t = omp_get_wtime();
    __lu_decomposition(N);
    printf("%lf\n", omp_get_wtime() - t);
    // _print_sq(l, N, 2);
    return 0;
}
