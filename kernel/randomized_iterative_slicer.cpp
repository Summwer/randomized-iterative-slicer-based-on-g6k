#include "randomized_iterative_slicer.h"


void save_matrix(ZZ_mat<mpz_t> mat, const char* filename){
  std::fstream myFile;
  std::cout<<filename<<std::endl;
  myFile.open(filename, ios_base::in|ios_base::out|ios_base::trunc);


  if (myFile.is_open()){
    // myFile << n << endl;
    // myFile << m << endl;
    myFile << mat <<endl;
    std::cout <<"Save file successful!" <<endl;
   
    myFile.close();
  }

}



//gh^2
FT compute_gh(std::vector<double> rr, unsigned int dim){
  // Compute the Gaussian Heuristic of the current block
  double const log_ball_square_vol = dim * std::log(M_PI) - 2.0 * std::lgamma(dim / 2.0 + 1);
  double log_lattice_square_vol = 0;
  for (unsigned int i = 0; i < dim; ++i)
  {
    log_lattice_square_vol += std::log(rr[i]);
  }
  FT gh = std::exp((log_lattice_square_vol - log_ball_square_vol) / (1.0 * dim));
  return gh;
}


Siever::Vec<VEC_ZT> Siever::update_vector(Vec<VEC_ZT> vec,unsigned int dim){
  vec.c = sim_hashes.compress(vec.v);
  vec.norm = std::inner_product(vec.v.begin(), vec.v.begin()+dim, vec.v.begin(),  static_cast<FT>(0.0));
  if(vec.norm < 0) {
    std::cout<<vec.norm;
    throw "The data type is too small to present the data!";
  }
  return vec;
}


Siever::Vec<VEC_ZT> Siever::recover_vector(std::array<ZT,MAX_SIEVING_DIM> x, Matrix<Z_NR<mpz_t>> B,int q, unsigned int dim){
  Vec<VEC_ZT> vec;
  for(unsigned int j = 0; j < dim; j++){
    vec.v[j] = 0; //Initiate
    for(unsigned int k = 0; k < dim; k++){
      vec.v[j] += x[k]*B[k][j].get_si();   
      // std::cout<<x[k]<<",";
      // std::cout<<B[k][j].get_si()<<",";
      // std::cout<<vec.v[j]<<" ";
    }
    vec.v[j] = mod_min(vec.v[j],q);
    // std::cout<<std::endl;
    
    //std::cout<<vec.v[j]<<" ";
  }  
  
  
  vec.close_x = x;
  vec = update_vector(vec,dim);
  return vec;
}
 
// Siever::Vec<mpz_t> Siever::recover_vector(std::array<ZT,MAX_SIEVING_DIM> x, Matrix<Z_NR<mpz_t>> B,int q, unsigned int dim){
//   Vec<mpz_t> vec;
  
//   for(unsigned int j = 0; j < dim; j++){
//     mpz_set_ui(vec.v[j] , 0); //Initiate
//     for(unsigned int k = 0; k < dim; k++){
//       mpz_t tmp1;
//       mpz_t tmp2;
//       B[k][j].get_mpz(tmp2);
//       mpz_mul_ui(tmp1,tmp2, x[k]);
//       mpz_add (vec.v[j], vec.v[j] , tmp1);   
//     }
//     vec.v[j] = mod_min(vec.v[j],q);
//   }  
//   vec.close_x = x;
//   vec = update_vector(vec,dim);
//   return vec;
// }






//the db after sieve is the coordinate of basis, thus we should recover them.
void Siever::recover_db(Matrix<Z_NR<mpz_t>> B, int q,unsigned int dim){
  vecs.resize(db.size());
  for(unsigned int j = 0; j < db.size(); j++){
    Entry e = db[cdb[j].i];
    vecs[j] = recover_vector(e.x, B, q, dim);
  } 
}


//the db after sieve is the coordinate of basis, thus we should recover them, recover an input db
void Siever::recover_db(std::vector<std::vector<ZT>> input_db, fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<FT>> M, int q,unsigned int dim){
  initialize_local_params(M,0,0, dim,dim);
  //std::cout<<M.b[1][0]<<std::endl;
  unsigned long N = input_db.size();
  reserve(N);
  cdb.resize(N);
  db.resize(N);
  vecs.resize(N);
  //std::array<ZT,MAX_SIEVING_DIM> x {};
  Entry e;
  for(unsigned int j = 0; j < input_db.size(); j++){
    for(unsigned int i = 0; i < dim; i++){
      e.x[i] = input_db[j][i];
    }  

    // if(j==db.size()-1){
    //   std::cout<<"x = ";
    //   for(unsigned int k = 0; k < dim; k++){
    //     for(unsigned int t = 0; t < dim; t++){
    //       //std::cout<<e.x[t]<<",";
    //       std::cout<<t<<","<<k<<",";
    //       std::cout<<M.b[t][k].get_si()<<" ";
    //     //std::cout<<vecs[db.size()-1].v[j]<<" ";
    //     }
    //     std::cout<<"--------------------"<<std::endl;
    //     break;
    //   }
      
    //   std::cout<<std::endl;
    // }
    
   
    //std::cout<<vecs[0].norm<<std::endl;
    //return;
    int large = 0;
    size_t np = dim / 2;
    np = (large > 3) ? 0 : np;
    //recompute_data_for_entry_babai<Recompute::recompute_all_and_consider_otf_lift>(e,np);
    recompute_data_for_entry<Recompute::recompute_all_and_consider_otf_lift>(e);

    //if (!uid_hash_table.insert_uid(e.uid)) continue;
    histo[histo_index(e.len)] ++;
    db[j] = e;
            
    CompressedEntry ce;
    ce.len = e.len;
    ce.c = e.c;
    ce.i = j;
    cdb[j] = ce;
    // std::cout<<ce.c<<std::endl;

    vecs[j] = recover_vector(e.x, M.b, q, dim);

    // if(j==db.size()-1){
    //   std::cout<<"x = ";
    //   for(unsigned int k = 0; k < dim; k++)
    //     std::cout<<e.x[k]<<" ";
    //   std::cout<<std::endl;
    // }
    

    
  } 
}


//sample  a vector t' in L(A+t), just sample a vector in L(A), and plus it to t.
Siever::Vec<VEC_ZT> Siever::sample_t_(Vec<VEC_ZT> target_vector, Matrix<Z_NR<mpz_t>> B, int q, unsigned int dim){
  Entry e;
  e = sample(0); //large = 0 
  //for(unsigned int i = 0; i < dim; i++){
  //  std::cout<<e.x[i]<<" ";
  //}
  
  Vec<VEC_ZT> t_ = recover_vector(e.x,B,q,dim);
  for(unsigned int i = 0; i < dim; i++){
    t_.v[i] = mod_min(target_vector.v[i] - t_.v[i],q);
  }//Contruct the t_ = L(A)+t 

  //recompute t_ in simhash and its norm. represent t as an Entry & CompressedEntry
  t_ = update_vector(t_,dim);
 
  return t_;
}


void Siever::initialize_local_params(fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<FT>> M,unsigned int ll, unsigned int l, unsigned int r, unsigned int dim){
    std::cout<<"--Initialize Local Parameters..."<<std::endl;
    //std::cout<<M.b[1][0]<<std::endl;
    FP_NR<FT> mu;
    FP_NR<FT> rri;
    std::vector<double> rr (dim,0);
    double* mu_double = new double[dim*dim]; //initialize the 1-dim array for mu
   
    for(unsigned int i=ll; i< r; i++){
        for(unsigned int j=ll; j< r; j++){
            //std::cout<<M.get_mu(mu,i,j)<<" ";
            if(i==j){
              rr[i] = M.get_r(rri,i,i).get_d();
              mu_double[i * dim + i] = rr[i];
            }else if(j<i){
              M.get_mu(mu,i,j);
              mu_double[i * dim + j] = mu.get_d();
            }
            else{
              mu_double[i * dim + j] = 0.0;
            }
        }
    }
    load_gso(dim, mu_double);
    initialize_local(ll,l,r);
    //std::cout<<M.b[1][0]<<std::endl;
    //3.initialize db
    //siever.shrink_db(0);
    reset_stats(); //clear all previous statistics
}


void Siever::call_gauss_sieve(fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<FT>> M,unsigned int ll, unsigned int l, unsigned int r, int q, unsigned int dim){
    std::cout<<"--Initialize Local Parameters..."<<std::endl;
    FP_NR<FT> mu;
    FP_NR<FT> rri;
    std::vector<double> rr (dim,0);
    double* mu_double = new double[dim*dim]; //initialize the 1-dim array for mu
   
    for(unsigned int i=ll; i< r; i++){
        for(unsigned int j=ll; j< r; j++){
            //std::cout<<M.get_mu(mu,i,j)<<" ";
            if(i==j){
              rr[i] = M.get_r(rri,i,i).get_d();
              mu_double[i * dim + i] = rr[i];
            }else if(j<i){
              M.get_mu(mu,i,j);
              mu_double[i * dim + j] = mu.get_d();
            }
            else{
              mu_double[i * dim + j] = 0.0;
            }
        }
    }
    load_gso(dim, mu_double);
    initialize_local(ll,l,r);
    
    //3.initialize db
    //siever.shrink_db(0);
    reset_stats(); //clear all previous statistics
  
    //4.sample vectors to grow db_size.
    //double db_size_factor = 3.2;
    //double db_size_base = pow(4/3.,1/2.);
    //int N = 500 + 10*dim + 2 * db_size_factor * pow(db_size_base, dim);//969;// pow(4/3.,dim/2);
    int N =  pow(4/3.,dim/2);
    //int large = 0;
    //siever.grow_db(ceil(N), large);
   
    //5.Start Gauss Sieve.
    std::cout<<"--Start Gauss Sieve..."<<std::endl;
    gauss_sieve(N); //如果提前找到足够数量的短向量，会提前跳出筛法循环，不一定会找到N多个短向量。
    std::cout<<"--Finish Gauss Sieve."<<std::endl;

    vecs.resize(db.size());
    recover_db(M.b,q,dim);
}

Siever::Vec<VEC_ZT> Siever::slicer(Vec<VEC_ZT> target_vector, Matrix<Z_NR<mpz_t>> B, FT norm_bound, int max_sample_times,int q, unsigned int dim){
  
  Siever::Vec<VEC_ZT> t_new1;
  Siever::Vec<VEC_ZT> t_new2;
  Siever::Vec<VEC_ZT> mint_ = target_vector;

  int sample_times = 1;
start_sample_t_:
  //3. sample  a vector t' in L(A+t), just sample a vector in L(A), and plus it to t. 
  Siever::Vec<VEC_ZT> t_ = sample_t_(target_vector, B,q,dim);
  
  // for(unsigned int i = 0; i < dim; i++){
  //   std::cout<< t_.v[i] << " ";
  // }
  // std::cout<<std::endl;
  // return target_vector;
  //4. sort db in order from nearby with target_vector to farway with target_vector. 
  //reduce the sample_vector by slicer algorithm, and get a short enough vector.
start_over:
  //std::cout<<t_.norm<<std::endl;
  for(unsigned int i = 0; i < db.size(); i++){
    if( t_.norm <= norm_bound){
        break;
    }
    
    if(UNLIKELY(is_reducible_maybe<XPC_THRESHOLD>(t_.c,vecs[i].c))){
        
      for(unsigned int j = 0; j < dim; j++){
        t_new1.v[j] = mod_min(t_.v[j] + vecs[i].v[j],q);
        t_new1.close_x[j] = mod_min(t_.close_x[j] - vecs[i].close_x[j],q);
      }
      t_new1 = update_vector(t_new1,dim);
  
      for(unsigned int  j = 0; j < dim; j++){
        t_new2.v[j] = mod_min(t_.v[j] - vecs[i].v[j],q);
        t_new2.close_x[j] = mod_min(t_.close_x[j] + vecs[i].close_x[j],q);
      }
      t_new2 = update_vector(t_new2,dim);

      if(t_new1.norm < t_.norm){
        t_ = t_new1;
        goto start_over;
      }
      else if(t_new2.norm < t_.norm){
        t_ = t_new2;
        goto start_over;
      }    
    }
  }

  if ( t_.norm < mint_.norm)
    mint_ = t_;
  //std::cout<<(sample_times+1)<<","<<mint_.norm<<std::endl;
  if( t_.norm > norm_bound and sample_times < max_sample_times){
    sample_times += 1;
    goto start_sample_t_;
  }

  //Siever::Vec<int> close_vector = recover_vector(mint_.close_x,B,q,dim);//t_ = t + close_vector
  
  // for(unsigned int i = 0; i < dim; i++){
  //     close_vector.v[i] = - close_vector.v[i];
  //     close_vector.close_x[i] = - close_vector.close_x[i];
  //     target_vector.close_x[i] = - close_vector.close_x[i];
  // }//t_ = t-close_vector

  // close_vector = update_vector(close_vector);
  if(sample_times > 1)
    std::cout<<"Sample Times = "<<sample_times<<std::endl;
  return mint_;

  // return close_vector;
}



vector<int> generate_random_poly(int q,unsigned int n){
  vector<int> poly (n,0);
  for(unsigned int i = 0; i< n; i++){
    //std::cout<<rand()<<" ";
    poly[i] = mod_min(rand() , q);
  }
  
  return poly;
}


//the poly contains the coefficients of m +/-1
vector<int> generate_bin_poly(unsigned int n,unsigned int m){
  vector<int> poly (n,0);
  for(unsigned int i = 0; i< m; i++){
    poly[i] = (rand() % 2) * 2 - 1 ;
  }

  std::random_shuffle(poly.begin(), poly.end());
  // for(unsigned int  i = 0; i< n; i++){
  //   std::cout<<poly[i]<<" ";
  // }
  return poly;
}


//poly a = a0 +a1x + ... + a(n-1)x^(n-1) --> [ a(n-2) a(n-3) ...  a0 -a(n-1)]
vector<int> blacklash(vector<int> poly){
  int element = - poly.front();
  poly.erase(poly.begin(),poly.begin()+1);
  poly.push_back(element);
  return poly;
}

//poly a = a0 +a1x + ... + a(n-1)x^(n-1)
//--> A = [[ a(n-1) a(n-2) ...  a1   a0   ]
//         [ a(n-2) a(n-3) ...  a0 -a(n-1)]
//         ....
//         [ a0    -a(n-1) ... -a2 -a1   ]]
Matrix<int> transform_poly_to_mat(vector<int> poly,unsigned int n){
  Matrix<int> Mat(n,n);
  reverse(poly.begin(),poly.end()); //[a0,a1,...a(n-1)] --> [a(n-1),a(n-2),...a0]
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < n; j++){
      Mat[i][j] = poly[j];
      //std::cout<<Mat[i][j]<<" ";
    }
    //std::cout<<std::endl;
    poly = blacklash(poly);
  }
  return Mat;
}

//return the abs minimum number
int mod_min(int a, int q){
  a = a % q;
  if(2 * a > q){
    a -= ceil(a/float(q))*q;
  } else if( 2 * a < -q){
    a -= floor(a/float(q))*q;
  }
  //std::cout<<a<<std::endl;
  return a;
}

// //return the abs minimum number
// mpz_t mod_min(mpz_t a, mpz_t q){
//   a = a % q;
//   if(2 * a > q){
//     a -= ceil(a/float(q))*q;
//   } else if( 2 * a < -q){
//     a -= floor(a/float(q))*q;
//   }
//   //std::cout<<a<<std::endl;
//   return a;
// }


//Calculate the result of poly = c*t
vector<int> mul(vector<int> c,vector<int> t, unsigned int n, int q){
  Matrix<int> C = transform_poly_to_mat(c,n);
  vector<int> ct (n,0);
  //std::cout<<std::endl;
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < n; j++){
      ct[i] += C[i][j]*t[j];
      //std::cout<<ct[i]<<" ";
    }
    // std::cout<<std::endl;
    // std::cout<<ct[i]<<std::endl;
    ct[i] = mod_min(ct[i], q);
    // std::cout<<ct[i]<<std::endl;
    // std::cout<<q<<std::endl;
  }
  reverse(ct.begin(),ct.end());
  return ct;
}


//lattice basis: [[qI 0 0]
//                [A -I 0]
//                [B 0 -I]] 
ZZ_mat<mpz_t> construct_MSIS_cvp_lattice_instance(int q,unsigned int n, unsigned int m){
  ZZ_mat<mpz_t> lattice_basis;
  std::vector<int> a;
  Matrix<int> A;
  
  lattice_basis.gen_zero((m+1)*n,(m+1)*n);
  for(unsigned int i = 0; i< (m+1)*n; i++){
    if(i<n)
      lattice_basis[i][i] = q;
    else
      lattice_basis[i][i] = -1;
  }

  
  for(unsigned int  iter = 0; iter < m; iter++){
    a = generate_random_poly(q,n); //poly = [a0,a1,...a(n-1)]
    A = transform_poly_to_mat(a,n); //poly --> blscklash matrix

    for(unsigned int i = 0; i < n; i++){
      for(unsigned int j = 0; j < n; j++){
        lattice_basis[i+(iter+1)*n][j] = A[i][j];
      }
    }
  }
  return lattice_basis;
}



//lattice basis: [[qI 0  0  0]
//                [A -I  0  0]
//                [B  0 -I  0]
//                [D  0  0 -I]] 
ZZ_mat<mpz_t> construct_MSIS_svp_lattice_instance(unsigned int n, unsigned int m, ZZ_mat<mpz_t> cvp_lattice,std::vector<int> ct){
  ZZ_mat<mpz_t> lattice_basis;
  lattice_basis.gen_zero((m+2)*n,(m+2)*n);
  Matrix<int> D = transform_poly_to_mat(ct,n);

  for(unsigned int i = 0; i < (m+2)*n; i++){
    for(unsigned int j = 0; j < (m+2)*n; j++){
      if(i < (m+1)*n and j < (m+1)*n)
        lattice_basis[i][j] = cvp_lattice[i][j];
      else{
        if( i >= (m+1)*n and j<n )
          lattice_basis[i][j] = D[i-(m+1)*n][j];
        else if(i==j)
          lattice_basis[i][j] = -1;
      }
    } 
  }
  return lattice_basis;
}




// int main(){
//     //test:Input lattice basis: Matrix A; Target vector: target_vector
//     // ZZ_mat<mpz_t> A;
//     // int status = read_file(A, "svpchallenge/svpchallengedim40seed0.txt");

//     int q = 4096;
//     unsigned int n = 20;
//     ZZ_mat<mpz_t> lattice_basis = construct_MSIS_cvp_lattice_instance(q,n);
    
//     //2.update gso and implent an LLL reduction
//     ZZ_mat<mpz_t> U;
//     ZZ_mat<mpz_t> UT;
//     MatGSO<Z_NR<mpz_t>, FP_NR<FT>> M(lattice_basis, U, UT, GSO_DEFAULT);
//     LLLReduction<Z_NR<mpz_t>, FP_NR<FT>> lll(M,LLL_DEF_DELTA, LLL_DEF_ETA, LLL_DEFAULT);
//     lll.lll(0,0,dim);
//     M.update_gso();

//     FT gamma = 1.05;
//     std::array<double,dim> rr;
//     FP_NR<FT> rri;
  
//     for(unsigned int i=0; i< dim; i++){
//       rr[i] = M.get_r(rri,i,i).get_d();
//     }
  
//     FT gh = compute_gh(rr);
//     FT target_norm = gamma * gh;
//     std::vector<int> c = generate_bin_poly(n,floor(n/2)); //the poly contains the coefficients of m +/-1
    
//     std::vector<int> c0 = generate_bin_poly(n,2); //the poly contains the coefficients of m +/-1
//     std::vector<int> t = generate_random_poly(q,n); //poly = [a0,a1,...a(n-1)]
  
//     std::vector<int> ct = mul(c,t,n);
//     std::vector<int> c0t = mul(c0,t,n);
    
//     std::cout<<"c0t = [ ";
//     for(unsigned int i = 0; i < n; i++){
//       std::cout<< c0t[i] << " ";
//     }
//     std::cout<<"]."<<std::endl;
    
    
//     std::cout<<"--Initialize class Siever--"<<std::endl;
//     //initiate siever
//     SieverParams params;
//     params.simhash_codes_basedir = "spherical_coding"; //directory of simhash spherical_coding files 
//     //std::cout<<params.threads<<std::endl;
//     unsigned long int seed = 0;
//     Siever siever(params,seed);

//     //2. Generate db for basis A
//     clock_t Gauss_siever_T0 = clock();
//     siever.call_gauss_sieve(M,0,0,M.d);
//     clock_t Gauss_siever_T = clock() - Gauss_siever_T0;
//     //store the recovered db.
//     siever.vecs.resize(siever.db.size());
 
//     siever.recover_db(M.b);

//     std::cout<<"--Finish Gauss Sieve, time cost = "<< Gauss_siever_T  * 1.0 / CLOCKS_PER_SEC <<"s--"<<std::endl;
//     std::cout<<siever.vecs[siever.db.size()-1].norm<<std::endl;


//     //3. Implement the randomized iterative slicer algorithm
//     Siever::Vec target_vector;
//     for(unsigned int i = 0; i< dim; i++){
//       if(i<n)
//         target_vector.v[i] = ct[i];
//       else
//         target_vector.v[i] = 0;
//       target_vector.close_x[i] = 0;
//     }
//     target_vector = siever.update_vector(target_vector);
//     std::cout<<"The origininal norm of target_vector is "<<target_vector.norm<<", "<<"the target norm is "<<target_norm<<"."<<std::endl;
//     //bool negative = true; //true: t_ = t + close_vector, the -close_vector is close to t;
//                           //false: t_ = t - close_vector, the close_vector is close to t.

//     std::cout<<"--Start slicer--"<<std::endl;
//     clock_t slicer_T0 = clock();
//     Siever::Vec t_ = siever.slicer(target_vector, M.b, target_norm,1000000); // return t_ = t - close_vector.
//     clock_t slicer_T = clock() - slicer_T0;

//     std::cout<<"The norm of reduced target_vector is "<<t_.norm<<"."<<std::endl;
//     std::cout<<"The reduced target_vector t_ = [ ";
//     for(unsigned int i = 0; i < dim; i++){
//       std::cout<< t_.v[i] << " ";
//     }
//     std::cout<<"]."<<std::endl;

//     std::cout<<"--Finish slicer, time cost = "<< slicer_T * 1.0 / CLOCKS_PER_SEC <<"s--"<<std::endl;
    
    
//     return 0;

// }