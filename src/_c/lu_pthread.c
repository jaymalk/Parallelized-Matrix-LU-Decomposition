#include "commons.h"
#include "lu_pthread.h"

int no_of_threads = 1;
double **m, **l, **u, *p;

typedef struct init_arg_struct {
    double ** mat;
    int size;
}init_arg_struct;

struct loop_arg {
    int start;
    int end;
    int extra1;
    int extra2;
};

struct init_loop_arg{
    double **mat;
    int start;
    int end;
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
void *swap_l(void *arguments);
void *lu(void *arguments);
void *mlu(void *arguments);
void __lu_decomposition(int size) {

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

        pthread_t thread_id[no_of_threads];
        struct loop_arg args;
        args.extra1 = k;
        args.extra2 = kf;
        
        int start_loop = 0;
        int end_loop = k;
        int gap = k/no_of_threads;
        for(int i=0;i<no_of_threads;i++){
            args.start = start_loop + i*gap;
            if(i==no_of_threads-1)
                args.end = end_loop;
            else
                args.end = start_loop+gap;

            pthread_create( &thread_id[i], NULL, &swap_l, (void *)&args);
        }

        for(int j=0; j < no_of_threads; j++)
            pthread_join( thread_id[j], NULL);
        
        // setting values
        u[k][k] = m[k][k];

        args.extra1 = k;
        start_loop = k+1;
        end_loop = size;
        gap = (size-k-1)/no_of_threads;
        for(int i=0;i<no_of_threads;i++){
            args.start = start_loop + i*gap;
            if(i==no_of_threads-1)
                args.end = end_loop;
            else
                args.end = start_loop+gap;

            pthread_create( &thread_id[i], NULL, &lu, (void *)&args);
        }

        for(int j=0; j < no_of_threads; j++)
            pthread_join( thread_id[j], NULL);
        
        args.extra1 = k;
        args.extra2 = size;
        start_loop = k+1;
        end_loop = size;
        gap = (size-k-1)/no_of_threads;
        for(int i=0;i<no_of_threads;i++){
            args.start = start_loop + i*gap;
            if(i==no_of_threads-1)
                args.end = end_loop;
            else
                args.end = start_loop+gap;

            pthread_create( &thread_id[i], NULL, &mlu, (void *)&args);
        }

        for(int j=0; j < no_of_threads; j++)
            pthread_join( thread_id[j], NULL);
    }
}


/*
 * 2D Matrix (Square), parallel initialisations.
 * @param (double ***): reference to 2D array
 * @param (int): order of the matrix
 */
void *threaded_alloc(void *arguments);
void *__init_2d(void *arguments) {
    struct init_arg_struct *args = (struct init_arg_struct*)arguments;
    double ** m = args->mat;
    int size = args->size;
    (m) = (double **)malloc(sizeof(double *)*size);

    pthread_t thread_id[no_of_threads];

    int start_loop = 0;
    int end_loop = size;
    int gap = size/no_of_threads;
    for(int i=0;i<no_of_threads;i++){
        struct init_loop_arg *arg = (struct init_loop_arg *)malloc(sizeof(struct init_loop_arg));
        arg->size = size;
        arg->mat = m;
        arg->start = start_loop + i*gap;
        if(i==no_of_threads-1)
            arg->end = end_loop;
        else
            arg->end = arg->start+gap;

        pthread_create( &thread_id[i], NULL, &threaded_alloc, (void *)arg);
    }

    for(int j=0; j < no_of_threads; j++)
        pthread_join( thread_id[j], NULL);
}


/*
 * Complete initialisations of matrices involved (parallely).
 * @param (double ***): matrix references
 * @param (int): order of matrices
 */

void *threaded_init(void *arguments);
void init(int N) {
    // 2D init
    pthread_t init_thread_id[3];
    {
        struct init_arg_struct *args = (struct init_arg_struct *)malloc(sizeof(struct init_arg_struct));
        args->size = N;
        args->mat = m;
        pthread_create( &init_thread_id[0], NULL, &__init_2d, (void *)args);
        printf("m %ld\n",sizeof(args->mat));
        args->mat = l;
        pthread_create( &init_thread_id[1], NULL, &__init_2d, (void *)args);
        
        args->mat = u;
        pthread_create( &init_thread_id[2], NULL, &__init_2d, (void *)args);
        
    }

    for(int j=0; j < 3; j++)
        pthread_join( init_thread_id[j], NULL);
 
    // 1D init
    p = (double *)malloc(sizeof(double)*N);

    pthread_t thread_id[no_of_threads];

    int start_loop = 0;
    int end_loop = N;
    int gap = N/no_of_threads;
    printf("Beginning random init %d %d %d\n",start_loop,end_loop,gap);
    for(int i=0;i<no_of_threads;i++){ 
        struct loop_arg *args = (struct loop_arg *)malloc(sizeof(struct loop_arg));
        args->extra1 = N;
    
        args->start = start_loop + i*gap;
        if(i==no_of_threads-1)
            args->end = end_loop;
        else
            args->end = args->start+gap;

        printf("Creating thread %d %d %d %d\n",i,args->start,args->end,args->extra1);
        pthread_create( &thread_id[i], NULL, &threaded_init, (void *)args);
    }

    for(int j=0; j < no_of_threads; j++)
        pthread_join( thread_id[j], NULL);
    printf("Ending Random Init\n");
    
}



int main(int argc, char const *argv[])
{
    int N = atoi(argv[1]);
    no_of_threads = atoi(argv[2]);
    
    // double t = omp_get_wtime();
    printf("Beginning Init\n");
    init(N);
    printf("Ending Init\n");
    // printf("Initialization %lf\n", omp_get_wtime() - t);
    
    // t = omp_get_wtime();
    printf("Beginning LU\n");
    __lu_decomposition(N);
    printf("Ending Init\n");
    // printf("%lf\n", omp_get_wtime() - t);
    // _print_sq(l, N, 2);
    return 0;
}


void *swap_l(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int k = args->extra1;
    int kf = args->extra2;
    for(int i=start; i<end; i++)
            swap_d(l[k]+i, l[kf]+i);
}

void *lu(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int k = args->extra1;
    for(int i=start; i<end; i++) {
        l[i][k] = m[i][k]/u[k][k];
        u[k][i] = m[k][i];
    }
}

void *mlu(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int k = args->extra1;
    int size = args->extra2;
    for(int i=start; i<end; i++) {
        for(int j=k+1; j<size; j++)
            m[i][j] -= (l[i][k]*u[k][j]);
    }
}

void *threaded_init(void *arguments){
    printf("Entering threaded_init ");
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int N = args->extra1;
    printf("Threaded_Init %d %d %d\n",start,end,N);
    for(int i=start; i<end; i++) {
        // perm.matrix
        printf("Updating P\n");
        p[i] = i;
        printf("Updated P %d\n",i);
        for(int j=0; j<N; j++) {
            // matrix
            printf("Updating m %d %d\n",i,j);
            m[i][j] = drand48();
            printf("Updated m\n");
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
    printf("Exiting threaded_init()\n");
}

void *threaded_alloc(void *arguments){
    struct init_loop_arg *args = (struct init_loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int size = args->size;
    double **m = args->mat;
    for(int i=start; i<end; i++){
        m[i] = (double *)malloc(sizeof(double)*size);
        printf("threaded alloc %d %ld\n",i,sizeof(m[i]));
    }
}