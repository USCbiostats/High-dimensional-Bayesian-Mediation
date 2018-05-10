/**
 * mcmc_mixture.c
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

/** Multi-variable Normal distribution */
#include <mvn.hpp>


using namespace hmlp;

int main( int argc, char *argv[] )
{
  using T = double;

  size_t d = 3;

  Data<T> mu( d, 1, 0 );
  Data<T> Sigma( d, d ); Sigma.randspd(); 

  MultiVariableNornal<T> my_mvn( mu, Sigma );

  return 0;
};
