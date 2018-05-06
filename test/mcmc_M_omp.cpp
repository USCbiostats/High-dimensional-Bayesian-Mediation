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
#include <bslmm_M.hpp>

#ifdef HMLP_MIC_AVX512
#include <hbwmalloc.h>
#endif

#define GFLOPS 1073741824 
#define TOLERANCE 1E-13

using namespace std;
using namespace hmlp;

int main( int argc, char *argv[] )
{
  using T = double;

  size_t n = 800;
  size_t w1 = 3;
  size_t w2 = 7;
  size_t q = 755;
  size_t q1 = 100;
  size_t q2 = 100;
  size_t burnIn = 30000;
  size_t niter = 50000;
  size_t permute = 1;

  if ( argc == 19 )
  {
    /** read parameters */
    sscanf( argv[ 1 ], "%lu", &n );
    sscanf( argv[ 2 ], "%lu", &w1 );
    sscanf( argv[ 3 ], "%lu", &w2 );
    sscanf( argv[ 4 ], "%lu", &q );
    sscanf( argv[ 5 ], "%lu", &burnIn );
    sscanf( argv[ 6 ], "%lu", &niter );
    sscanf( argv[ 7 ], "%lu", &q1 );
    sscanf( argv[ 8 ], "%lu", &q2 );
    sscanf( argv[ 9 ], "%lu", &permute );
  }
  else
  {
    printf( "\n[usage] ./mcmc.x <n> <w> <q> <q1> <niter>\n\n" );
  }


  string Y_filename(       argv[ 10 ] );
  string M_filename(       argv[ 11 ] );
  string A_filename(       argv[ 12 ] );
  string C1_filename(      argv[ 13 ] );
  string C2_filename(      argv[ 14 ] ); 
  string beta_m_filename(  argv[ 15 ] );
  string alpha_a_filename( argv[ 16 ] );
  string pi_m_filename(    argv[ 17 ] );
  string pi_a_filename(    argv[ 18 ] );

  #pragma omp parallel for 
  for ( int p = 1; p < permute; p ++ )
  {
    Data<T> Y( n, 1 );
    Data<T> M( n, q ); 
    Data<T> A( n, 1 );
    Data<T> C1( n, w1 );
    Data<T> C2( n, w2 );
    Data<T> beta_m( 1, q );
    Data<T> alpha_a( 1, q ); 
    Data<T> pi_m( 1, q );
    Data<T> pi_a( 1, q );

    string Y_name     = Y_filename                                        + string( ".txt" );
    string M_name     = M_filename       + string( "_" ) + to_string( p ) + string( ".txt" );
    string A_name     = A_filename                                        + string( ".txt" );
    string C1_name    = C1_filename                                       + string( ".txt" );
    string C2_name    = C2_filename                                       + string( ".txt" );
    string beta_name  = beta_m_filename  + string( "_" ) + to_string( p ) + string( ".txt" );
    string alpha_name = alpha_a_filename + string( "_" ) + to_string( p ) + string( ".txt" );
    string pi_m_name  = pi_m_filename    + string( "_" ) + to_string( p ) + string( ".txt" );
    string pi_a_name  = pi_a_filename    + string( "_" ) + to_string( p ) + string( ".txt" );

    Y.readmatrix(        n,  1,     Y_name );
    M.readmatrix(        n,  q,     M_name );
    A.readmatrix(        n,  1,     A_name );
    C1.readmatrix(       n, w1,    C1_name );
    C2.readmatrix(       n, w2,    C2_name );
    beta_m.readmatrix(   1,  q,  beta_name );
    alpha_a.readmatrix(  1,  q, alpha_name );
    pi_m.readmatrix(     1,  q,  pi_m_name );
    pi_a.readmatrix(     1,  q,  pi_a_name );

    mcmc::mcmc<T>( Y, A, M, C1, C2, beta_m, alpha_a, pi_m, pi_a, 
        n, w1, w2, q, q1, q2, burnIn, niter, p );
  }

  return 0;
};
