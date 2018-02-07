/**
 * mcmc.c
 *
 * Chenhan D. Yu
 *
 * Department of Computer Science, University of Texas at Austin
 *
 * Purpose: 
 *
 * Todo:
 *
 * Modification:
 *
 **/

#include <tuple>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include <hmlp_blas_lapack.h>

#include <data.hpp>
#include <bslmm.hpp>

#ifdef HMLP_MIC_AVX512
#include <hbwmalloc.h>
#endif

#define GFLOPS 1073741824 
#define TOLERANCE 1E-13

using namespace hmlp::mcmc;

int main( int argc, char *argv[] )
{
  using T = double;

  size_t n = 1230;
  size_t w1 = 3;
  size_t w2 = 13;
  size_t q = 7661;
  size_t q1 = 100;
  size_t q2 = 100;
  size_t burnIn = 30000;
  size_t niter = 50000;

	if ( argc == 16 )
	{
		/** read parameters */
		sscanf( argv[ 1 ], "%lu", &n );
		sscanf( argv[ 2 ], "%lu", &w1 );
		sscanf( argv[ 3 ], "%lu", &w2 );
		sscanf( argv[ 4 ], "%lu", &q );
		sscanf( argv[ 5 ], "%lu", &burnIn );
		sscanf( argv[ 6 ], "%lu", &niter );
	}
	else
	{
		printf( "\n[usage] ./mcmc.x <n> <w> <q> <q1> <niter>\n\n" );
	}


    std::string Y_filename(       argv[  7 ] );
    std::string M_filename(       argv[  8 ] );
    std::string A_filename(       argv[ 9 ] );
    std::string C1_filename(      argv[ 10 ]  );
	std::string C2_filename(      argv[ 11 ] ); 
    std::string beta_m_filename(  argv[ 12 ]  );
    std::string alpha_a_filename( argv[ 13 ] );
    std::string pi_m_filename(    argv[ 14 ]  );
    std::string pi_a_filename(    argv[ 15 ]  );

	//std::string Y_filename( "bmi3.txt" );
    //std::string M_filename( "sig_shore.txt" );
    //std::string A_filename( "race_st.txt" );
    //std::string C1_filename( "age.sex.txt" );
	//std::string C2_filename( "age.sex.10pc.txt" );

	//std::string beta_m_filename( "beta_m2.txt" );
    //std::string alpha_a_filename( "alpha_a2.txt" );
    //std::string pi_m_filename( "pi_m2.txt" );
    //std::string pi_a_filename( "pi_a2.txt" );

    hmlp::Data<T> Y( n, 1 );
    hmlp::Data<T> M( n, q ); 
    hmlp::Data<T> A( n, 1 );
    hmlp::Data<T> C1( n, w1 );
	hmlp::Data<T> C2( n, w2 );

    hmlp::Data<T> beta_m( 1, q );
    hmlp::Data<T> alpha_a( 1, q ); 
    hmlp::Data<T> pi_m( 1, q );
    hmlp::Data<T> pi_a( 1, q );

    Y.readmatrix( n, 1, Y_filename );
    M.readmatrix( n, q, M_filename );
    A.readmatrix( n, 1, A_filename );
    C1.readmatrix( n, w1, C1_filename );
	C2.readmatrix( n, w2, C2_filename );

	beta_m.readmatrix( 1, q, beta_m_filename );
    alpha_a.readmatrix( 1, q, alpha_a_filename );
    pi_m.readmatrix( 1, q, pi_m_filename );
    pi_a.readmatrix( 1, q, pi_a_filename );

	//hmlp::Data<T> X( 2, 3 ); X.randn();
	//X.WriteFile( "X.m" );

  mcmc<T>( Y, A, M, C1, C2, beta_m, alpha_a, pi_m, pi_a, n, w1, w2, q, q1, q2, burnIn, niter );

  return 0;
};
