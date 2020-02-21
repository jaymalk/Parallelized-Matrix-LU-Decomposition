
#include "include.h"

#include "commons.h"


/*
 * Swapping double pointed values.
 * @param _1, _2 (double *)
 */
void swap_d(double * _1, double * _2) {
    double _t = *_1;
    *_1 = *_2; 
    *_2 = _t;
}

/*
 * Swapping int pointed values.
 * @param _1, _2 (int *)
 */
void swap_i(int * _1, int * _2) {
    int _t = *_1;
    *_1 = *_2; 
    *_2 = _t;
}


/*
 * Swapping double pointed rows.
 * @param _1, _2 (double **)
 */
void swap_d_r(double ** _1, double ** _2) {
    double* _t = *_1;
    *_1 = *_2; 
    *_2 = _t;
}

/*
 * Standard Matrix Multiplication.
 * @param _1, _2(double **): matrices involved
 * @param result(double **): result stored in this matrix
 * @param size(int): Size of the square matrices involved
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
void __matmul(double ** _1, double ** _2, double ** result, int size) {
    printf("Multiplying Matrix\n");
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++){           
            double temp = 0.0;
            for(int k=0; k<size; k++){
                temp += (_1[i][k] * _2[k][j]);
            }
            result[i][j] = temp;
        }
    printf("Matrix multiplication done\n");
}

/*
 * Standard Matrix Multiplication.
 * @param original, result(double **): matrices involved
 * @param result(double **): result stored in this matrix
 * @param size(int): Size of the square matrices involved
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
double checker(double **original, double ** result, int *p, int size, int _fd){
    double error = 0.0;
    int r;
    for(int i=0; i<size; i++){
        r=p[i];
        for(int j=0; j<size; j++)
            error += (original[r][j] - result[i][j])*(original[r][j] - result[i][j]);
    }
    return error;
}


/*
 * Reading a matrix from a file.
 * @param filename (const char *): name of the file (local path)
 * @param matrix (double **): container matrix
 * @param size (int): matrix order
 */
void read_matrix(const char * filename, double ** matrix, int size) {
    // Opening the file
    FILE *_file = fopen(filename, "r+");
    // Checking for error
    if(_file == NULL) {
        fprintf(__stderrp, "Can't read the matrix, from the file, %s", filename);
        exit(2);
    }
    // Reading the matrix array
    for(uint16_t i=0; i<size; i++)
        for(uint16_t j=0; j<size; j++)
           if ( fscanf(_file, "%lf", &(matrix[i][j])));
    // Closing the file
    fclose(_file);
}


/*
 * Writing a matrix to a file.
 * @param filename (const char *): name of the file (local path)
 * @param matrix (double **): container matrix
 * @param size (int): matrix order
 */
void write_matrix(const char * filename, double ** matrix, int size) {
    // Opening the file
    FILE *_file = fopen(filename, "w+");
    // Checking for error
    if(_file == NULL) {
        fprintf(__stderrp, "Can't write the matrix, to the file, %s", filename);
        exit(2);
    }
    // Reading the matrix array
    for(uint16_t i=0; i<size; i++) {
        for(uint16_t j=0; j<size; j++)
            fprintf(_file, "%6.2lf", matrix[i][j]);
        fprintf(_file, "\n");
    }
    // Closing the file
    fclose(_file);
}


/*
 * Write a vector to a file.
 * @param filename (const char *): name of the file (local path)
 * @param matrix (int *): container vector
 * @param size (int): matrix order
 */
void write_vector(const char * filename,  int * vector, int size) {
    // Opening the file
    FILE *_file = fopen(filename, "w+");
    // Checking for error
    if(_file == NULL) {
        fprintf(__stderrp, "Can't write the vector, to the file, %s", filename);
        exit(2);
    }
    // Reading the matrix array
    for(uint16_t i=0; i<size; i++)
        fprintf(_file, "%d\n", vector[i]);
    // Closing the file
    fclose(_file);
}