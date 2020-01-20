
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
 * Swapping double pointed rows.
 * @param _1, _2 (double **)
 */
void swap_d_r(double ** _1, double ** _2) {
    double* _t = *_1;
    *_1 = *_2; 
    *_2 = _t;
}


/*
 * Ordered printing of a square-matrix.
 * @param _mat (double **): matrix
 * @param _sze (int): size of the matrix
 * @param _fd (file directory no.)
 */
void _print_sq(double **_mat, int _sze, int _fd) {
    char *s = (char *)malloc(10);
    for(int i=0; i<_sze; i++) {
        for(int j=0; j<_sze; j++) {
            sprintf(s, "%16.12lf", _mat[i][j]);
            write(_fd, s, 10);
        }
        write(_fd, "\n", 1);
    }
}