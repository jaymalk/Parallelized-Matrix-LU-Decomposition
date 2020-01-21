
#include "include.hpp"
#include "commons.hpp"
#include "lu_openmp.hpp"

using namespace std;

int no_of_threads = 1;
/*
 * The LU Decomposition Function. (LOWER)
 * @param a_ (MAT): input array (COPY)
 * @param l_, u_ (MAT): output lu decomposition
 * @param p_ (vector<int>): output permutation matrix
 */
void __lu_decomposition(MAT a_, MAT& l_, MAT& u_, std::vector<int>& p_) {

    // Local Variables
    double max; int kf, size = p_.size();

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
            cerr << "Singular Matrix...\nExiting\n";
            exit(2);
        }

        // swapping
        swap(p_[k], p_[kf]);
        swap(a_[k], a_[kf]);

#       pragma omp parallel for num_threads(no_of_threads)   
        for(int i=0; i<k; i++)
            swap(l_[k][i], l_[kf][i]);
        
        // setting values
        u_[k][k] = a_[k][k];
#       pragma omp parallel for num_threads(no_of_threads)
        for(int i=k+1; i<size; i++) {
            l_[i][k] = a_[i][k]/u_[k][k];
            u_[k][i] = a_[k][i];
        }

#       pragma omp parallel for num_threads(no_of_threads)
        for(int i=k+1; i<size; i++)
#           pragma omp parallel for
            for(int j=k+1; j<size; j++)
                a_[i][j] -= (l_[i][k]*u_[k][j]);
    }
}

/*
 * Complete initialisations of matrices involved.
 * @param (MAT): matrices
 * @param (int): size
 */
void init(MAT& m_, MAT& l_, MAT& u_, std::vector<int>& p_, int _size) {
    
    // Initialisations
#   pragma omp parallel sections num_threads(4)
    {
#       pragma omp section
        m_ = vector<vector<double>>(_size, vector<double>(_size, 0));
#       pragma omp section
        u_ = vector<vector<double>>(_size, vector<double>(_size, 0));
#       pragma omp section
        l_ = vector<vector<double>>(_size, vector<double>(_size, 0));
#       pragma omp section
        p_ = vector<int>(_size, 0);
    }

    // Fillings
#   pragma omp parallel for num_threads(no_of_threads)
    for(int i=0; i<_size; i++) {
        // permutation
        p_[i] = i;
#       pragma omp parallel for
        for(int j=0; j<_size; j++) {
            // matrix
            m_[i][j] = _rand(_seed);
            // u & l (conditional)
            if(i == j) {
                l_[i][j] = 1;
                u_[i][j] = 1;
            }
            else if (i > j)     l_[i][j] = 1;
            else                u_[i][j] = 1;
        }
    }
}



int main(int argc, char const *argv[])
{
    int N = stoi(argv[1]);
    no_of_threads = stoi(argv[2]);
    MAT m, l, u;
    vector<int> p;
    double t = omp_get_wtime();
    init(m, l, u, p, N);
    printf("Initialization %lf\n", omp_get_wtime() - t);
    t = omp_get_wtime();
    // __print(m);
    __lu_decomposition(m, l, u, p);
    printf("%lf\n", omp_get_wtime() - t);
    cout << "done\n";
    // __print(l);
    // cout << endl;
    // __print(u);
    // cout << endl;
    // __print_permute(m, p);
    // cout << endl;
    // MAT x = __matmul(l, u);
    // __print(x);
    return 0;
}
