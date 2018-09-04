
#include <tuple>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <limits>

/** Multi-variable Normal distribution */
#include <mvn.hpp>


using namespace hmlp;

int main( int argc, char *argv[] )
{
  using T = double;

  size_t d = 3;
  size_t n_mixtures = 4;

  vector<Data<T>> Sigma_Mixture( n_mixtures );
  for ( size_t i = 0; i < n_mixtures; i ++ ) 
    Sigma_Mixture[ i ].resize( d, d );

  Data<T> Sigma0 = Sigma_Mixture[ 0 ];
  Data<T> Sigma1 = Sigma_Mixture[ 1 ];


  Data<T> mu( d, 1, 0 );
  Data<T> Sigma( d, d ); Sigma.randspd<false>(); 

  MultiVariableNormal<T> my_mvn( mu, Sigma );

  Sigma.Print();
  auto det = my_mvn.Determinant();
  auto logdet = my_mvn.LogDeterminant();
  auto sample = my_mvn.SampleFrom( 10 );
  sample.Print();


  return 0;
};
