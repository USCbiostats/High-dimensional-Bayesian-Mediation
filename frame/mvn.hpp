#include <hmlp_blas_lapack.h>
#include <data.hpp>

using namespace std;

namespace hmlp
{

template<typename T>
class MultiVariateNormal
{
  public:

    MultiVariateNormal( Data<T> mu, Data<T> &Sigma )
    {
      this->d = mu.row();
      this->mu = mu;
      this->Sigma = Sigma;
      /** Cholesky factorization (POTRF): Sigma = LL' */
      xportrf( "Lower", d, Sigma.data(), d );
      /** Compute the determinant from the Cholesky factorization */
      for ( uint64_t i = 0; i < d; i ++ ) det *= Sigma( i, i ) * Sigma( i, i );
    };

    T Determinant() { return ret; };

    T LogDeterminant() { return std::log( ret ) }

    Data<T> SampleFrom( uint64_t num_of_samples )
    {
      Data<T> X( d, num_of_samples );
      /** TODO: Normal( 0, 1 ) */
      X.randn();
      /** Compute L * X using TRMM. */
      xtrmm( "Left", "Lower", "No Transpose", "Not Unit", d, num_of_samples,
          1.0, Sigma.data(), d, X.data(), d )
      return X;
    };

  private:

    /** Dimension (number of variables) */
    uint64_t d = 0;

    /** d-by-1 expectation */
    Data<T> mu;

    /** d-by-d variance-covariance matrix */
    Data<T> Sigma;

    T det = 1;

};


}; /** end namespace hmlp */
