#include <iostream>
#include <stack>
#include "siever.h"
#include "../g6k-fplll/tests/test_utils.h"
#include "../g6k-fplll/fplll/gso_interface.h"
#include "../g6k-fplll/fplll/lll.h"
#include "../g6k-fplll/fplll/defs.h"
#include <time.h>


#ifndef DIM
#define DIM 80
#endif

FT compute_gh(std::vector<double> rr, unsigned int dim);
void save_matrix(ZZ_mat<mpz_t> mat, const char* filename);
vector<int> generate_random_poly(int q,unsigned int n); //Generate a random poly in Zq/(X^n+1)

//the poly contains the coefficients of m +/-1 in Zq/(X^n+1)
vector<int> generate_bin_poly(unsigned int n,unsigned int m);

//poly a = a0 +a1x + ... + a(n-1)x^(n-1) --> [ a(n-2) a(n-3) ...  a0 -a(n-1)]
vector<int> blacklash(vector<int> poly);

//poly a = a0 +a1x + ... + a(n-1)x^(n-1)
//--> A = [[ a(n-1) a(n-2) ...  a1   a0   ]
//         [ a(n-2) a(n-3) ...  a0 -a(n-1)]
//         ....
//         [ a0    -a(n-1) ... -a2 -a1   ]]
Matrix<int> transform_poly_to_mat(vector<int> poly,unsigned int n);

//Calculate the result of poly = c*t
vector<int> mul(vector<int> c,vector<int> t, unsigned int n, int q);
int mod_min(int a, int q);
//lattice basis: [[qI 0 0]
//                [A -I 0]
//                [B 0 -I]] 
ZZ_mat<mpz_t> construct_MSIS_cvp_lattice_instance(int q,unsigned int n,unsigned int m);

//lattice basis: [[qI 0  0  0]
//                [A -I  0  0]
//                [B  0 -I  0]
//                [D  0  0 -I]] 
ZZ_mat<mpz_t> construct_MSIS_svp_lattice_instance(unsigned int n, unsigned int m, ZZ_mat<mpz_t> cvp_lattice,std::vector<int> ct);