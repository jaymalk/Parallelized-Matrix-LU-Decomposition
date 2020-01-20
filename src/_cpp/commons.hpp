
#include "include.hpp"

#ifndef __COM_PP
#define __COM_PP


/*
 * Ordered printing of a square-matrix, on stderr.
 * @param _mat (vector<vector<double>>): matrix
 */
void __print(std::vector<std::vector<double>>&);


/*
 * Ordered printing of a square-matrix, keeping in mind the permutation.
 * @param _mat (vector<vector<double>>): matrix
 */
void __print_permute(std::vector<std::vector<double>>& _mat, std::vector<int>& _p);

/*
 * Standard Matrix Multiplication.
 * @param (std::vector<std::vector<double>>&, std::vector<std::vector<double>>&): matrices involved
 * 
 * NOTE: Both matrices are assumed to be square and of the same order.
 */
std::vector<std::vector<double>> __matmul(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);


#endif /* __COM */