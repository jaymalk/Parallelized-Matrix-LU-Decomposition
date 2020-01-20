
#include "include.hpp"

#ifndef __LU_OMP_PP
#define __LU_OMP_PP

#define MAT std::vector<std::vector<double>>

// Random Generator
std::random_device __rd;
std::mt19937 _seed(__rd());
std::uniform_real_distribution<> _rand(0.0, 1.0);


/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (double **): input array (COPY)
 * @param l_, u_ (double **): output lu decomposition
 * @param p_ (double *): output permutation matrix
 * 
 * NOTE: Copy of matrix (a_) must be passed as it gets overwritten.
 */
void __lu_decomposition(MAT, MAT&, MAT&, std::vector<int>&);


/*
 * Complete initialisations of matrices involved.
 * @param (MAT): matrices
 * @param (int): size
 */
void init(MAT&, MAT&, MAT&, std::vector<int>&, int);



#endif /* __LU_OMP_PP */