#include "commons.h"
#include "lu_pthread.h"

int no_of_threads = 1;
double **m, **mcopy, **l, **u;
int *p;

// Struct to initialize matrices
typedef struct init_arg_struct {
    double *** mat;
    int size;
}init_arg_struct;

//struct to parallelise loops
struct loop_arg {
    int start;
    int end;
    int extra1;
    int extra2;
};

//Struct to initialize matrices
struct init_loop_arg{
    double ***mat;
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
        swap_i(p+k, p+kf);
        swap_d_r(m+k, m+kf);

        //Parallelsing loop to swap values of l[k] and l[k']
        pthread_t swap_l_thread_id[no_of_threads];
        
        int start_loop = 0;
        int end_loop = k;
        int gap = k/no_of_threads;
        for(int i=0;i<no_of_threads;i++){
            //A new struct for every part of loop
            struct loop_arg *args = (struct loop_arg *)malloc(sizeof(struct loop_arg));
            args->extra1 = k;
            args->extra2 = kf;
            args->start = start_loop + i*gap;
            if(i==no_of_threads-1)
                args->end = end_loop;
            else
                args->end = args->start+gap;

            pthread_create( &swap_l_thread_id[i], NULL, &swap_l, (void *)args);
        }

        for(int j=0; j < no_of_threads; j++)
            pthread_join( swap_l_thread_id[j], NULL);
        
        // setting values
        u[k][k] = m[k][k];

        //Parallelsing loop to compute values of l and u
        pthread_t lu_thread_id[no_of_threads];

        start_loop = k+1;
        end_loop = size;
        gap = (size-k-1)/no_of_threads;
        for(int i=0;i<no_of_threads;i++){
            //A new struct for every part of loop
            struct loop_arg *args = (struct loop_arg *)malloc(sizeof(struct loop_arg));
            args->extra1 = k;
            args->start = start_loop + i*gap;
            if(i==no_of_threads-1)
                args->end = end_loop;
            else
                args->end = args->start+gap;
            
            pthread_create( &lu_thread_id[i], NULL, &lu, (void *)args);
        }

        //Joining Threads
        for(int j=0; j < no_of_threads; j++)
            pthread_join( lu_thread_id[j], NULL);
        

        //Parallelsing loop to swap values of a
        pthread_t mlu_thread_id[no_of_threads];
        start_loop = k+1;
        end_loop = size;
        gap = (size-k-1)/no_of_threads;
        for(int i=0;i<no_of_threads;i++){
            //A new struct for every part of loop
            struct loop_arg *args = (struct loop_arg *)malloc(sizeof(struct loop_arg));
            args->extra1 = k;
            args->extra2 = size;
            args->start = start_loop + i*gap;
            if(i==no_of_threads-1)
                args->end = end_loop;
            else
                args->end = args->start+gap;

            pthread_create( &mlu_thread_id[i], NULL, &mlu, (void *)args);
        }

        //Joining Threads
        for(int j=0; j < no_of_threads; j++)
            pthread_join( mlu_thread_id[j], NULL);
    }
}


/*
 * 2D Matrix (Square), parallel initialisations.
 * @param (double ***): reference to 2D array
 * @param (int): order of the matrix
 */
void *threaded_alloc(void *arguments);
void *__init_2d(void *arguments) {

    //Extracting parameters from struct
    struct init_arg_struct *args = (struct init_arg_struct*)arguments;
    double *** _m = args->mat;
    int size = args->size;
    (*_m) = (double **)malloc(sizeof(double *) * size);

    //Parallelsing loop to allocate space to each row of matrix
    pthread_t thread_id[no_of_threads];

    int start_loop = 0;
    int end_loop = size;
    int gap = size/no_of_threads;
    for(int i=0;i<no_of_threads;i++){
        //A new struct for every part of loop
        struct init_loop_arg *arg = (struct init_loop_arg *)malloc(sizeof(struct init_loop_arg));
        arg->size = size;
        arg->mat = _m;
        arg->start = start_loop + i*gap;
        if(i==no_of_threads-1)
            arg->end = end_loop;
        else
            arg->end = arg->start+gap;

        pthread_create( &thread_id[i], NULL, &threaded_alloc, (void *)arg);
    }

    //Joining loops
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
    pthread_t* init_thread_id;
    init_thread_id = malloc(4*sizeof(pthread_t));
    {
        //Allocating memory to m
        struct init_arg_struct *args_m = (struct init_arg_struct *)malloc(sizeof(struct init_arg_struct));
        args_m->size = N;
        args_m->mat = &m;
        pthread_create( &init_thread_id[0], NULL, &__init_2d, (void *)args_m);
        
        //Allocating memory to mcopy to make copy of matrix a
        struct init_arg_struct *args_mcopy = (struct init_arg_struct *)malloc(sizeof(struct init_arg_struct));
        args_mcopy->size = N;
        args_mcopy->mat = &mcopy;
        pthread_create( &init_thread_id[1], NULL, &__init_2d, (void *)args_mcopy);
        
        //Allocating memory to l
        struct init_arg_struct *args_l = (struct init_arg_struct *)malloc(sizeof(struct init_arg_struct));
        args_l->size = N;
        args_l->mat = &l;
        pthread_create( &init_thread_id[2], NULL, &__init_2d, (void *)args_l);

        //Allocating memory to u       
        struct init_arg_struct *args_u = (struct init_arg_struct *)malloc(sizeof(struct init_arg_struct));
        args_u->size = N;
        args_u->mat = &u;
        pthread_create( &init_thread_id[3], NULL, &__init_2d, (void *)args_u);
    }
    //Joining Threads
    for(int j=0; j < 4; j++)
        pthread_join( init_thread_id[j], NULL);

    p = (int *)malloc(sizeof(int) * N);

    //Parallelising loop to randomly initialize the matrices
    pthread_t thread_id[no_of_threads];

    int start_loop = 0;
    int end_loop = N;
    int gap = N/no_of_threads;
    for(int i=0;i<no_of_threads;i++){ 
        struct loop_arg *args = (struct loop_arg *)malloc(sizeof(struct loop_arg));
        args->extra1 = N;
    
        args->start = start_loop + i*gap;
        if(i==no_of_threads-1)
            args->end = end_loop;
        else
            args->end = args->start+gap;

        pthread_create( &thread_id[i], NULL, &threaded_init, (void *)args);
    }

    //Joining Threads
    for(int j=0; j < no_of_threads; j++)
        pthread_join( thread_id[j], NULL);
    
}



int main(int argc, char const *argv[])
{
    //Extract parameters
    int N = atoi(argv[1]);
    no_of_threads = atoi(argv[2]);
    
    //Computing time for Initializing the matrices
    double t = omp_get_wtime();
    init(N);
    printf("Time taken for Initialization %lf\n", omp_get_wtime() - t);

    // Reading the matrices from the file
    read_matrix(argv[3], m, N);
    read_matrix(argv[3], mcopy, N);
    
    //Start counting time for LU decomposition, m now redundant
    t = omp_get_wtime();
    __lu_decomposition(N);
    printf("Time taken for LU Decomposition %lf\n", omp_get_wtime() - t);

    // Matrix multiplication
    __matmul(l,u,m,N);
    // Calculating norm
    printf("The L(2,1) norm is: %lf \n", checker(mcopy, m, p, N, 2));
    // Writing l and u
    write_matrix("L.txt", l, N);
    write_matrix("U.txt", u, N);
    return 0;
}


//Pointer function for pthread to swap values of l[k] and l[k']
void *swap_l(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int k = args->extra1;
    int kf = args->extra2;
    for(int i=start; i<end; i++)
            swap_d(l[k]+i, l[kf]+i);
}

//Pointer function for pthread to compute values of l and u
void *lu(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int k = args->extra1;
    for(int i=start; i<end; i++) {
        l[i][k] = (m)[i][k] / (u)[k][k];
        u[k][i] = m[k][i];
    }
}

//Pointer function for pthread to compute values of matrix a
void *mlu(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int k = args->extra1;
    int size = args->extra2;
    for(int i=start; i<end; i++) {
        for(int j=k+1; j<size; j++)
            m[i][j] -= (l[i][k] * u[k][j]);
    }
}

//Pointer function for pthread to randomly initialize the matrices
void *threaded_init(void *arguments){
    struct loop_arg *args = (struct loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int N = args->extra1;
    for(int i=start; i<end; i++) {
        p[i] = i;
        for(int j=0; j<N; j++) {
            m[i][j] = drand48();
            mcopy[i][j] = m[i][j];
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

//Pointer function for pthread to allocate memory to each row
void *threaded_alloc(void *arguments){
    struct init_loop_arg *args = (struct init_loop_arg*)arguments;
    int start = args->start;
    int end = args->end;
    int size = args->size;
    double ***m = args->mat;
    for(int i=start; i<end; i++){
        (*m)[i] = (double *)malloc(sizeof(double)*size);
    }
}