
#include "include.hpp"
#include "commons.hpp"


/*
 * Ordered printing of a square-matrix.
 * @param _mat (vector<vector<double>>): matrix
 */
void __print(std::vector<std::vector<double>>& _mat) {
    int _sze = _mat.size();
    for(int i=0; i<_sze; i++) {
        for(int j=0; j<_sze; j++) {
            std::cerr << std::setw(15) << std::setprecision(8) << _mat[i][j];
        }
        std::cerr << "\n";
    }
}

/*
 * Ordered printing of a square-matrix, keeping in mind the permutation.
 * @param _mat (vector<vector<double>>): matrix
 */
void __print_permute(std::vector<std::vector<double>>& _mat, std::vector<int>& _p) {
    int sze = _mat.size(), r;
    for(int i=0; i<sze; i++) {
        r = _p[i];
        for(int j=0; j<sze; j++) {
            std::cerr << std::setw(15) << std::setprecision(8) << _mat[r][j];
        }
        std::cerr << "\n";
    }
}

/*
 * Standard Matrix Multiplication.
 * @param (std::vector<std::vector<double>>&, std::vector<std::vector<double>>&): matrices involved
 * 
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
std::vector<std::vector<double>> __matmul(std::vector<std::vector<double>>& _1, std::vector<std::vector<double>>& _2) {
    int size = _1.size();
    std::vector<std::vector<double>> _m = std::vector<std::vector<double>>(size, std::vector<double>(size, 0));
    for(int i=0; i<size; i++)
        for(int j=0; j<size; j++)
            for(int k=0; k<size; k++)
                _m[i][j] += (_1[i][k] * _2[k][j]);
    return _m;
}