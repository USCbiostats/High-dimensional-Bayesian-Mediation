#ifndef MCMC_HPP
#define MCMC_HPP

#include <assert.h>
#include <typeinfo>
#include <algorithm>
#include <random>
#include <limits>
#include <cstddef>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;
using namespace hmlp;


namespace hmlp
{
namespace mcmc
{

template <typename RealType = double>
class beta_distribution
{
  public:
    typedef RealType result_type;

    class param_type
    {
      public:
        typedef beta_distribution distribution_type;

        explicit param_type(RealType a = 2.0, RealType b = 2.0)
          : a_param(a), b_param(b) { }

        RealType a() const { return a_param; }
        RealType b() const { return b_param; }

        bool operator==(const param_type& other) const
        {
          return (a_param == other.a_param &&
              b_param == other.b_param);
        }

        bool operator!=(const param_type& other) const
        {
          return !(*this == other);
        }

      private:
        RealType a_param, b_param;
    };

    explicit beta_distribution(RealType a = 2.0, RealType b = 2.0)
      : a_gamma(a), b_gamma(b) { }
    explicit beta_distribution(const param_type& param)
      : a_gamma(param.a()), b_gamma(param.b()) { }

    void reset() { }

    param_type param() const
    {
      return param_type(a(), b());
    }

    void param(const param_type& param)
    {
      a_gamma = gamma_dist_type(param.a());
      b_gamma = gamma_dist_type(param.b());
    }

    template <typename URNG>
    result_type operator()(URNG& engine)
    {
      return generate(engine, a_gamma, b_gamma);
    }

    template <typename URNG>
    result_type operator()(URNG& engine, const param_type& param)
    {
      gamma_dist_type a_param_gamma(param.a()),
                      b_param_gamma(param.b());
      return generate(engine, a_param_gamma, b_param_gamma); 
    }

    result_type min() const { return 0.0; }
    result_type max() const { return 1.0; }

    result_type a() const { return a_gamma.alpha(); }
    result_type b() const { return b_gamma.alpha(); }

    bool operator==(const beta_distribution<result_type>& other) const
    {
      return (param() == other.param() &&
          a_gamma == other.a_gamma &&
          b_gamma == other.b_gamma);
    }

    bool operator!=(const beta_distribution<result_type>& other) const
    {
      return !(*this == other);
    }

  private:
    typedef std::gamma_distribution<result_type> gamma_dist_type;

    gamma_dist_type a_gamma, b_gamma;

    template <typename URNG>
    result_type generate(URNG& engine,
        gamma_dist_type& x_gamma,
        gamma_dist_type& y_gamma)
    {
      result_type x = x_gamma(engine);
      return x / (x + y_gamma(engine));
    }
};

  template <typename CharT, typename RealType>
    std::basic_ostream<CharT>& operator<<(std::basic_ostream<CharT>& os,
        const beta_distribution<RealType>& beta)
    {
      os << "~Beta(" << beta.a() << "," << beta.b() << ")";
      return os;
    }

  template <typename CharT, typename RealType>
    std::basic_istream<CharT>& operator>>(std::basic_istream<CharT>& is,
        beta_distribution<RealType>& beta)
    {
      std::string str;
      RealType a, b;
      if (std::getline(is, str, '(') && str == "~Beta" &&
          is >> a && is.get() == ',' && is >> b && is.get() == ')') {
        beta = beta_distribution<RealType>(a, b);
      } else {
        is.setstate(std::ios::failbit);
      }
      return is;
    }



template<typename T>
void Normalize( hmlp::Data<T> &A )
{
  size_t m = A.dim();
  size_t n = A.num();

  for ( size_t j = 0; j < n; j ++ )
  {
    T mean = 0.0;
    T stde = 0.0;

    /** mean */
    for ( size_t i = 0; i < m; i ++ ) mean += A( i, j );
    mean /= m;

    /** standard deviation */
    for ( size_t i = 0; i < m; i ++ ) 
    {
      T tmp = A( i, j ) - mean;
      stde += tmp * tmp;
    }
    stde /= m;
    stde = std::sqrt( stde );

    /** normalize */
    for ( size_t i = 0; i < m; i ++ ) 
    {
      A[ j * m + i ] = ( A( i, j ) - mean ) / stde;
    }
  }

}; /** end Normalize() */


template<typename T>
class Variables
{
  public:

	Variables(
		hmlp::Data<T> &userY,
		hmlp::Data<T> &userA,
		hmlp::Data<T> &userM,
		hmlp::Data<T> &userC1,
		hmlp::Data<T> &userC2,
 	  hmlp::Data<T> &userbeta_m,
		hmlp::Data<T> &useralpha_a,
		hmlp::Data<T> &userpi_m,
		hmlp::Data<T> &userpi_a,
		size_t n, size_t w1, size_t w2, size_t q, size_t q1, size_t q2, size_t permute )
	  : Y( userY ), A( userA ), M( userM ), C1( userC1 ), C2( userC2 ),
	    beta_m( userbeta_m ), alpha_a( useralpha_a ), pi_m( userpi_m ), pi_a( userpi_a )

	{
    this->n = n;
    this->w1 = w1;
    this->w2 = w2;
    this->q = q;
    this->q1 = q1;
    this->q2 = q2;
    this->permute = permute;
    
    /** Initialize my_samples here. */
    my_samples.resize( 1000, 6 * q + 5, 0.0 );

    /** generate synthetic data Y, A, M and E */
    //beta_a.resize( 1, 1, 0.5 );
    //beta_m.resize( 1, q, 0.0 ); beta_m.randn( 0.0, 2.0 );
    //for ( size_t i = q1; i < beta_m.size(); i ++ ) beta_m[ i ] = 0.0;
    //A.resize( n, 1 ); A.randn( 0.0, 1.0 );
    //M.resize( n, q ); M.randn( 0.0, 2.0 );
    //E.resize( n, 1 ); E.randn( 0.0, 1.0 );

    /** Y <- beta_a*A + M %*% beta_m + E */
    //Y = E;
    //for ( size_t i = 0; i < n; i ++ ) Y[ i ] += beta_a[ 0 ] * A[ i ];
    //xgemm( "N", "N", n, 1, q, 1.0, M.data(), n, beta_m.data(), q, 1.0, Y.data(), n );

    /** normalize all synthetic data */
    //Normalize( Y );
    //Normalize( A );
    //Normalize( M );

    /** beta distribution */
    //beta_distribution<T> dist_pi_m( um, vm );
    //beta_distribution<T> dist_pi_a( ua, va );
    //pi_m.resize( 1, q, 0.0);
    //pi_a.resize( 1, q, 0.0);
    //for ( size_t i = 0; i < q; i ++) 
    //{
    // pi_m[ i ]  = dist_pi_m( generator );
    // pi_a[ i ]  = dist_pi_a( generator );
    //};

    //printf( "here\n" ); fflush( stdout );

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed( seed );

    /** gamma distribution initialization */
    std::gamma_distribution<T>  dist_m0(  km0, 1.0 / lm0 );     
    std::gamma_distribution<T>  dist_m1(  km1, 1.0 / lm1 );
    std::gamma_distribution<T>  dist_a(  ka, 1.0 / la );
    std::gamma_distribution<T> dist_ma0( kma0, 1.0 / lma0 ); 
    std::gamma_distribution<T> dist_ma1( kma1, 1.0 / lma1 );
    std::gamma_distribution<T>  dist_g(  kg, 1.0 / lg );
    std::gamma_distribution<T>  dist_e(  ke, 1.0 / le );
    sigma_m0  = 1.0 / dist_m0 ( generator ); 
    sigma_m1  = 1.0 / dist_m1 ( generator );
    sigma_a  = 1.0 / dist_a ( generator );
    sigma_ma0 = 1.0 / dist_ma0( generator );
    sigma_ma1 = 1.0 / dist_ma1( generator );
    sigma_g  = 1.0 / dist_g ( generator );
    sigma_e  = 1.0 / dist_e ( generator );

    /** assign initial values */
    r1.resize( 1, q, 0.0 );
    r3.resize( 1, q, 0.0 );

    /** resize beta_a */
    beta_a.resize( 0, 0 ); 
    beta_a.resize( 1, 1 ); 
    //beta_a.randn( 0, std::sqrt( sigma_a ) );
    beta_a[ 0 ] = 0.0135567741;

    //beta_m.resize( 0, 0 ); 
    //beta_m.resize( 1, q, 0.0 );
    //alpha_a.resize( 0, 0 ); 
    //alpha_a.resize( 1, q, 0.0 );

    beta_c.resize( 1, w1, 0.0 );
    alpha_c.resize( w2, q, 0.0 );

    /** compute column 2-norm */
    A2norm.resize( 1, 1, 0.0 );
    for ( size_t i = 0; i < n; i ++ ) A2norm[ 0 ] += A[ i ] * A[ i ];

    //printf( "A2norm %.3E \n", A2norm[ 0 ] ); fflush( stdout );
    
    M2norm.resize( 1, q, 0.0 );
    for ( size_t j = 0; j < q; j ++ )
      for ( size_t i = 0; i < n; i ++ )
        M2norm[ j ] += M[ j * n + i ] * M[ j * n + i ];

    //printf( "M2norm %.3E \n", M2norm[ 10 ] ); fflush( stdout );    

    C1_2norm.resize( 1, w1, 0.0);
    for ( size_t j = 0; j < w1; j ++ )
      for ( size_t i = 0; i < n; i ++ )
        C1_2norm[ j ] += C1[ j*n + i ] * C1[ j*n + i ];

    C2_2norm.resize( 1, w2, 0.0);
    for ( size_t j = 0; j < w2; j ++ )
      for ( size_t i = 0; i < n; i ++ )
        C2_2norm[ j ] += C2[ j*n + i ] * C2[ j*n + i ];			

  };

  void Residual( size_t it )
  {
    /** res1 = Y - M * beta_m - beta_a * A */
    res1 = Y;
    for ( size_t i = 0; i < n; i ++ )
      res1[ i ] -= beta_a[ 0 ] * A[ i ];
    xgemm( "N", "N", n, 1, q, -1.0, M.data(), n, 
        beta_m.data(), q, 1.0, res1.data(), n ); 

    //xgemm( "N", "N", n, 1, w1, -1.0, C1.data(), n,  
    //    beta_c.data(), w1, 1.0, res1.data(), n );

    /** res2 = M - A * alpha_a - C * alpha_c */
    res2 = M;
    xgemm( "N", "N", n, q, 1, -1.0, A.data(), n, 
        alpha_a.data(), 1, 1.0, res2.data(), n );

    //xgemm( "N", "N", n, q, w2, -1.0, C2.data(), n, 
    //    alpha_c.data(), w2, 1.0, res2.data(), n ); 

    res2_c = M;
    //xgemm( "N", "N", n, q, w2, -1.0, C2.data(), n,  
    //   alpha_c.data(), w2, 1.0, res2_c.data(), n );
  };

	 T PostDistribution1( hmlp::Data<T> &my_pi_m )
	 {
	     T logc1 = 0.0;
		   for ( size_t j = 0; j < q; j ++ )
		   {
		     logc1 += std::log ( my_pi_m[ j ] ) * r1[ j ] + std::log ( 1 - my_pi_m[ j ] ) * ( 1 - r1[ j ] );
		   }
       return logc1;
   };

   T PostDistribution2( hmlp::Data<T> &my_pi_a )
   {
     T logc2 = 0.0;
     for ( size_t j = 0; j < q; j ++ )
     {
       logc2 += std::log ( my_pi_a[ j ] ) * r3[ j ] + std::log ( 1 - my_pi_a[ j ] ) * ( 1 - r3[ j ] );
     }
     return logc2;
   };


   T old = 0.0;

   size_t count = 0;

   void Iteration( size_t burnIn, size_t it )
   { 

     if ( it % 30000 == 0 ) 
     { 
       printf( "Iter %4lu sigma_m0 %.3E sigma_e %.3E sigma_g %.3E sigma_ma0 %.3E sigma_m1 %.3E sigma_a %.3E sigma_ma1 %.3E \n", 
           it, sigma_m0, sigma_e, sigma_g, sigma_ma0, sigma_m1, sigma_a, sigma_ma1 ); fflush( stdout ); 
     }
     /** Update res1, res2, res2_c */
     if ( it == 0 ) Residual( it );
     //printf( "Iter %4lu res1 %.3E res2 %.3E \n", it, res1[ 0 ], res2[ 0 ] ); fflush( stdout ); 

     /** var_m = 1.0 / ( 1 / sigma_m + M2norm / sigma_e ) */
     var_m0.resize( 1, q, 0.0 );
     var_m1.resize( 1, q, 0.0 );
     for ( size_t j = 0; j < q; j ++ )
     {
       var_m0[ j ] = 1.0 / ( 1.0 / sigma_m0 + M2norm[ j ] / sigma_e );
       var_m1[ j ] = 1.0 / ( 1.0 / sigma_m1 + M2norm[ j ] / sigma_e );
     }
 
     //printf( "Iter %4lu var_m0 %.3E var_m1 %.3E \n", it, var_m0[ 0 ], var_m1[ 0 ] ); fflush( stdout ); 

     /** sigma_e and sigma_g */
     T ke1 = ke + n / 2.0;
     T kg1 = q * ( n / 2.0 ) + kg;
     T le1 = 0.0;
     T lg1 = 0.0;
     for ( size_t i = 0; i < n; i ++ )
       le1 += res1[ i ] * res1[ i ];
     for ( size_t i = 0; i < n * q; i ++ )
       lg1 += res2[ i ] * res2[ i ];
     le1 = 1.0 / ( le1 / 2.0 + le );
     lg1 = 1.0 / ( lg1 / 2.0 + lg );
     std::gamma_distribution<T> dist_e( ke1, le1 );
     std::gamma_distribution<T> dist_g( kg1, lg1 );
     sigma_e  = 1.0 / dist_e( generator );
     sigma_g  = 1.0 / dist_g( generator );

     /** var_alpha_a and var_a */
     var_alpha_a0.resize( 1, 1, 0.0 ); var_alpha_a1.resize( 1, 1, 0.0 );
     var_a.resize( 1, 1, 0.0 );
     var_alpha_a0[ 0 ] = sigma_g / ( sigma_g / sigma_ma0 + A2norm[ 0 ] );
     var_alpha_a1[ 0 ] = sigma_g / ( sigma_g / sigma_ma1 + A2norm[ 0 ] );
     var_a[ 0 ] = sigma_e / ( sigma_e / sigma_a + A2norm[ 0 ] ); 

     //printf( "Iter %4lu sigma_e %.3E var_a %.3E \n", it, sigma_e, var_a[ 0 ] ); fflush( stdout ); 

     for ( size_t j = 0; j < q; j ++ )
     {
       /** update res1, res2 */

       /** mu_mj, mu_alpha_aj */
       T mu_mj = 0.0;
       T mu_alpha_aj = 0.0;

       for ( size_t i = 0; i < n; i ++ )
       {
         mu_mj += M( i, j ) * ( res1[ i ] + M( i, j ) * beta_m[ j ] );
         mu_alpha_aj += A[ i ] * ( res2_c[ j * n + i ] );
       }
       T mu_mj0 = mu_mj / ( ( sigma_e / sigma_m0 ) + M2norm[ j ] ); 
       T mu_mj1 = mu_mj / ( ( sigma_e / sigma_m1 ) + M2norm[ j ] );
       T mu_alpha_aj0 = mu_alpha_aj * ( var_alpha_a0[ 0 ] / sigma_g );
       T mu_alpha_aj1 = mu_alpha_aj * ( var_alpha_a1[ 0 ] / sigma_g );
       
       //printf( "Iter %4lu q %4lu mu_mj0 %.3E mu_mj1 %.3E \n", it, j, mu_mj0, mu_mj1 ); fflush( stdout ); 
       /** beta_m[ j ] = r1[ j ] * randn( mu_mj, var_m[ j ] ) */
       old = beta_m[ j ];
       std::normal_distribution<T> dist_norm_m0( mu_mj0, std::sqrt( var_m0[ j ] ) );
       std::normal_distribution<T> dist_norm_m1( mu_mj1, std::sqrt( var_m1[ j ] ) );
       beta_m[ j ] = r1[ j ] * dist_norm_m1( generator ) + 
         ( 1 - r1[ j ]) *dist_norm_m0( generator );
       
       //printf( "Iter %4lu q %4lu r1 %.3E beta_m %.3E \n", it, j, r1[ j ], beta_m[ j ] ); fflush( stdout ); 

       for ( size_t i = 0; i < n; i ++ )
       {
         res1[ i ] = res1[ i ] + ( old - beta_m[ j ] ) * M( i, j );
       }

       //printf( "Iter %4lu q %4lu res1 %.3E res1 %.3E \n", it, j, res1[ 0 ], res1[ 1 ] ); fflush( stdout ); 

       /** alpha_a[ j ] = r3[ j ] * randn( mu_alpha_aj, var_alpha_a ) */
       old = alpha_a[ j ];
       std::normal_distribution<T> dist_alpha_a0( mu_alpha_aj0, std::sqrt( var_alpha_a0[ 0 ] ) ); 
       std::normal_distribution<T> dist_alpha_a1( mu_alpha_aj1, std::sqrt( var_alpha_a1[ 0 ] ) );
       alpha_a[ j ] = r3[ j ] * dist_alpha_a1( generator ) +
       ( 1 - r3[ j ] ) * dist_alpha_a0( generator );
       for ( size_t i = 0; i < n; i ++ )
       {
       res2[ j*n + i ] = res2[ j*n + i ] + ( old - alpha_a[ j ] ) * A[ i ];
       }

       /** update r1[ j ] */
       T const2 = mu_mj1 * mu_mj1 / ( 2 * var_m1[ j ] ) - mu_mj0 * mu_mj0 / ( 2 * var_m0[ j ] ) + 
         0.5 * std::log( var_m1[ j ] / sigma_m1 ) - 0.5 * std::log( var_m0[ j ] / sigma_m0 ) + 
         std::log( pi_m[ j ] / ( 1.0 - pi_m[ j ] ) );
       if ( const2 < 300 )
       {
         const2 = std::exp( const2 );
         //std::default_random_engine generator;
         std::bernoulli_distribution distribution( const2 / (1.0+const2));
         if ( distribution (generator) )
           r1[ j ] = 1.0;
         else
           r1[ j ] = 0.0;
       }
       else
       {
         r1[ j ] = 1.0;
       }

       //printf( "Iter %4lu q %4lu r1 %.3E pi_m %.3E \n", it, j, r1[ j ], pi_m[ j ] ); fflush( stdout ); 

       /** update r3[ j ] */
       T const5 = mu_alpha_aj1 * mu_alpha_aj1 / ( 2 * var_alpha_a1[ 0 ] ) - 
       mu_alpha_aj0 * mu_alpha_aj0 / ( 2 * var_alpha_a0[ 0 ] ) + 
       0.5 * std::log( var_alpha_a1[ 0 ] / sigma_ma1 ) - 0.5 * std::log( var_alpha_a0[ 0 ] / sigma_ma0 ) + std::log( pi_a[ j ] / ( 1.0 - pi_a[ j ] ) );
       if ( const5 < 300 )
       {
         const5 = std::exp( const5 );
         //std::default_random_engine generator;
         std::bernoulli_distribution distribution( const5 / (1.0+const5));
         if ( distribution (generator) )
           r3[ j ] = 1.0;
         else
           r3[ j ] = 0.0;
       }
       else
       {
         r3[ j ] = 1.0;
       }

       beta_distribution<T> dist_pi_m( um + r1[ j ], vm + 1 - r1[ j ] );
       beta_distribution<T> dist_pi_a( ua + r3[ j ], va + 1 - r3[ j ] ); 
       pi_m[ j ]  = dist_pi_m( generator );
       pi_a[ j ]  = dist_pi_a( generator );

       //printf( "Iter %4lu q %4lu r1 %.3E pi_m %.3E \n", it, j, r1[ j ], pi_m[ j ] ); fflush( stdout ); 
       //for ( size_t j1 = 0; j1 < w2; j1 ++ )
       // {
       //   T mu_alpha_cj = 0.0;
       //   old = alpha_c( j1, j );
       //   for ( size_t i = 0; i < n; i ++ )
       //   {
       //     mu_alpha_cj += C2( i, j1 ) * ( res2[ j*n + i ] + alpha_c( j1, j ) * C2( i, j1 )  );
       //   }
       //   mu_alpha_cj = mu_alpha_cj / ( C2_2norm[ j1 ] ); 
       //   std::normal_distribution<T> dist_alpha_c( mu_alpha_cj, std::sqrt( sigma_g / C2_2norm[ j1 ] ) );
       //   alpha_c( j1, j )  = dist_alpha_c( generator );
       //   for ( size_t i = 0; i < n; i ++ )
       //   {
       //     res2[ j*n + i ] = res2[ j*n + i ] + ( old - alpha_c( j1, j ) ) * C2( i, j1 );
       //     res2_c[ j*n + i ] = res2_c[ j*n + i ] + ( old - alpha_c( j1, j ) ) * C2( i, j1 );
       //   }

       //} /** end for each j1 < w2 */

     } /** end for each j < q */

     //for ( size_t j = 0; j < w1; j ++)
     //{
     //T mu_cj = 0.0;
     //old = beta_c[ j ];
     //for ( size_t i = 0; i < n; i ++ )
     //{
     //  mu_cj += C1( i, j ) * ( res1[ i ] + C1( i, j ) * beta_c[ j ]);
     //  mu_cj = mu_cj / C1_2norm[ j ];
     //  std::normal_distribution<T> dist_beta_c( mu_cj, std::sqrt( sigma_e / C1_2norm[ j ] ) );
     //    beta_c[ j ]  = dist_beta_c( generator );
     //}
     //for ( size_t i = 0; i < n; i ++ )
     //{
     //   res1[ i ] = res1[ i ] + ( old - beta_c[ j ] ) * C1( i, j );
     //}

     //}


     /** update beta_a */
     T mu_a = 0.0;
     old = beta_a[ 0 ];
     for ( size_t i = 0; i < n; i ++ )
       mu_a += A[ i ] * ( res1[ i ] + beta_a[ 0 ] * A[ i ] );
     mu_a *= ( var_a[ 0 ] / sigma_e );
     std::normal_distribution<T> dist_beta_a( mu_a, std::sqrt( var_a[ 0 ] ) );
     beta_a[ 0 ] = dist_beta_a( generator );
     for ( size_t i = 0; i < n; i ++ )
     {
       res1[ i ] = res1[ i ] + ( old - beta_a[ 0 ] ) * A[ i ];
     }

     //printf( "Iter %4lu res1 %.3E res1 %.3E \n", it, res1[ 0 ], res1[ 1 ] ); fflush( stdout ); 

     /** update sigma_m, sigma_a and sigma_ma */
     T const6 = 0.0;
     T const7 = 0.0;
     T const8 = 0.0;
     T const9 = 0.0;

     for ( size_t j = 0; j < q; j ++ )
      {
        const6 += r1[ j ];
        const7 += beta_m[ j ] * beta_m[ j ] * r1[ j ];
        const8 += r3[ j ];
        const9 += alpha_a[ j ] * alpha_a[ j ] * r3[ j ];
      }
      const6 /= 2.0;
      const7 /= 2.0;
      const8 /= 2.0;
      const9 /= 2.0;

      std::gamma_distribution<T>  dist_m1( const6 +  km1, 1.0 / ( const7 +  lm1 ) );
      std::gamma_distribution<T>  dist_a(    0.5 +  ka, 1.0 / ( beta_a[ 0 ] * beta_a[ 0 ] / 2.0 + la ) );
      std::gamma_distribution<T> dist_ma1( const8 + kma1, 1.0 / ( const9 + lma1 ) );
      sigma_m1  = 1.0 / dist_m1 ( generator );
      sigma_a  = 1.0 / dist_a ( generator );
      sigma_ma1 = 1.0 / dist_ma1 ( generator );

      //printf( "Iter %4lu sigma_a %.3E sigma_m1 %.3E \n", it, sigma_a, sigma_m1 ); fflush( stdout );       

      const6 = 0.0;
      const7 = 0.0;
      const8 = 0.0;
      const9 = 0.0;

	    for ( size_t j = 0; j < q; j ++ )
      {
        const6 += 1 - r1[ j ];
        const7 += beta_m[ j ] * beta_m[ j ] * ( 1 - r1[ j ] );
        const8 += 1 - r3[ j ];
        const9 += alpha_a[ j ] * alpha_a[ j ] * ( 1 - r3[ j ] );
      }
      const6 /= 2.0;
      const7 /= 2.0;
      const8 /= 2.0;
      const9 /= 2.0;

      std::gamma_distribution<T>  dist_m0( const6 +  km0, 1.0 / ( const7 +  lm0 ) );
      std::gamma_distribution<T> dist_ma0( const8 + kma0, 1.0 / ( const9 + lma0 ) );
      sigma_m0  = 1.0 / dist_m0 ( generator );
      sigma_ma0 = 1.0 / dist_ma0 ( generator );
      
      //printf( "Iter %4lu sigma_a %.3E sigma_m0 %.3E \n", it, sigma_a, sigma_m0 ); fflush( stdout ); 
     //hmlp::Data<T> log_my_pi_m( 1, q ); log_my_pi_m.rand( -0.01, 0.01 );
	   //hmlp::Data<T> log_my_pi_a( 1, q ); log_my_pi_a.rand( -0.01, 0.01 );
     //hmlp::Data<T> my_pi_m( 1, q );
     //hmlp::Data<T> my_pi_a( 1, q );

      //T probab = 0.0;
	    //hmlp::Data<T> my_unif( 1, 1 );
	    //for ( size_t j = 0; j < q; j ++ )
	    //{
		    //my_pi_m[ j ] = pi_m[ j ] * std::exp( log_my_pi_m[ j ] );
		    //my_pi_m[ j ] = std::abs( my_pi_m[ j ] );
		    //if ( my_pi_m[ j ] > 1.0 ) my_pi_m[ j ] = 1.0 / my_pi_m[ j ];
		    //if ( my_pi_m[ j ] < 1.0/q ) my_pi_m[ j ] = 1.0 / ( q*q*my_pi_m[ j ] );
	   
        //my_pi_a[ j ] = pi_a[ j ]* std::exp( log_my_pi_a[ j ] );
        //my_pi_a[ j ] = std::abs( my_pi_a[ j ] );
        //if ( my_pi_a[ j ] > 1.0 ) my_pi_a[ j ] = 1.0 / my_pi_a[ j ];
        //if ( my_pi_a[ j ] < 1.0/q ) my_pi_a[ j ] = 1.0 / ( q*q*my_pi_a[ j ] );
      //}

      //probab = PostDistribution1( my_pi_m ) -
      //         PostDistribution1(    pi_m );
      //my_unif.rand( 0.0, 1.0 );
      //if ( probab > std::log (my_unif[ 0 ]) )
      //{

        //for ( size_t j = 0; j < q; j ++ )
        //{	
          //pi_m[ j ] = my_pi_m[ j ];
        //}

      //}

      //probab = PostDistribution2( my_pi_a ) -
      //         PostDistribution2(    pi_a );
      //my_unif.rand( 0.0, 1.0 );
      //if ( probab > std::log (my_unif[ 0 ]) )
      //{

        //for ( size_t j = 0; j < q; j ++ )
        //{
          //pi_a[ j ] = my_pi_a[ j ];
        //}

      //}


      if ( it > burnIn && it % 50 == 0 )
      {
        //std::ofstream outfile;
        //std::string outfilename = std::string( "results_pY_" ) + std::to_string( (int)q1 ) + std::to_string( (int)q2 ) + std::string( "_" ) + std::to_string( (int)permute ) + std::string( ".txt" );
        //outfile.open( outfilename.data(), std::ios_base::app );

        //std::ofstream outfile;
        //outfile.open("results_shore.txt", std::ios_base::app);
        for ( size_t i = 0; i < q; i +=1 ) 
        {
          my_samples( count, 6*i   ) = beta_m[ i ];
          my_samples( count, 6*i+1 ) = pi_m[ i ];
          my_samples( count, 6*i+2 ) = (int)r1[ i ];
        
          my_samples( count, 6*i+3 ) = alpha_a[ i ];
          my_samples( count, 6*i+4 ) = pi_a[ i ];
          my_samples( count, 6*i+5 ) = (int)r3[ i ];

        }		  

        my_samples( count, 6*q ) = beta_a[ 0 ]; 
        my_samples( count, 6*q + 1 ) = sigma_m0;
        my_samples( count, 6*q + 2 ) = sigma_m1;
        my_samples( count, 6*q + 3 ) = sigma_ma0;
        my_samples( count, 6*q + 4 ) = sigma_ma1;
        count += 1;

      }

      if ( count >= 999 ) 
      {
        string my_samples_filename = string( "results_" ) + to_string( (int)q1 ) + to_string( (int)q2 ) + string( ".txt" );
        my_samples.WriteFile( my_samples_filename.data() );
      }

   };


    size_t n;

    size_t w1;

  	size_t w2;

    size_t q;

    size_t q1;

    size_t q2;

    size_t permute;

    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator;

    T km0  = 2.0;

    T lm0  = 1E-6;

  	T km1 = 2.0;

  	T lm1 = 1.0;

    T ka  = 2.0;

    T la  = 1.0;
  
    T kma0 = 2.0;

    T lma0 = 1E-6;

	  T kma1 = 2.0;

	  T lma1 = 1.0;
   
    T ke  = 2.0;

    T le  = 1.0;

    T kg  = 2.0;

    T lg  = 1.0;

    T um  = 1.0;

    T ua  = 1.0;

    T vm  = 2.0;

    T va  = 2.0;

    T sigma_m0, sigma_m1;

    T sigma_a;

    T sigma_ma0, sigma_ma1;

    T sigma_g;

    T sigma_e;

    hmlp::Data<T> &alpha_a;

    hmlp::Data<T> beta_a;

    hmlp::Data<T> &beta_m;

    hmlp::Data<T> &pi_m, &pi_a;

	  hmlp::Data<T> alpha_c;

	  hmlp::Data<T> beta_c;

    hmlp::Data<T> &A;

    hmlp::Data<T> A2norm;

    hmlp::Data<T> &M;

    hmlp::Data<T> M2norm;

    hmlp::Data<T> &Y;

	  hmlp::Data<T> &C1;

	  hmlp::Data<T> &C2;

	  hmlp::Data<T> C1_2norm;

	  hmlp::Data<T> C2_2norm;

    hmlp::Data<T> r1;

    hmlp::Data<T> r3;

    hmlp::Data<T> my_samples;
    /** in Iteration() */

    hmlp::Data<T> res1;
    hmlp::Data<T> res2;
    hmlp::Data<T> res2_c;

    hmlp::Data<T> var_m0, var_m1;

    hmlp::Data<T> var_alpha_a0, var_alpha_a1;

    hmlp::Data<T> var_a;

  private:
};

template<typename T>
void mcmc( hmlp::Data<T> &Y,
	         hmlp::Data<T> &A,
		       hmlp::Data<T> &M,
		       hmlp::Data<T> &C1,
		       hmlp::Data<T> &C2,
		       hmlp::Data<T> &beta_m,
	         hmlp::Data<T> &alpha_a,
		       hmlp::Data<T> &pi_m,
		       hmlp::Data<T> &pi_a,
	         size_t n, size_t w1, size_t w2, size_t q, size_t q1, size_t q2, size_t burnIn, size_t niter, size_t permute )
{
  Variables<T> variables( Y, A, M, C1, C2, beta_m, alpha_a, pi_m, pi_a, n, w1, w2, q, q1, q2, permute );

  std::srand(std::time(nullptr));
  
  for ( size_t it = 0; it < niter; it ++ )
  {
    variables.Iteration( burnIn, it );
  }

};


};
};

#endif // ifndef MCMC_HPP
