#include <VectorAdd.hpp>

using namespace std;

namespace helloworld
{

/** float version */
vector<float> VectorAdd_f32( vector<float> &a, vector<float> &b )
{
  return VectorAdd<float>( a, b );
};


/** double version */
vector<double> VectorAdd_f64( vector<double> &a, vector<double> &b )
{
  return VectorAdd<double>( a, b );
};

/** integer version */
vector<int> VectorAdd_d32( vector<int> &a, vector<int> &b )
{
  return VectorAdd<int>( a, b );
};

}; /** end namespace helloword */
