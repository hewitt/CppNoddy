/// \file TrivialComplex.cpp
/// \ingroup Tests
/// \ingroup Generic
/// A sanity test for std::complex class.

#include <Types.h>

namespace CppNoddy
{
  namespace Example
  {
    /// A std::complex function for z^3 - 1.
    inline D_complex Cfn( const D_complex& z )
    {
      return z * z * z - 1.;
    }
  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Complex: simple test ============================\n";
  cout << "\n";

  D_complex Z( 0.0, 0.0 );
  Z = Example::Cfn( D_complex( 0.5, 0.5 ) );

  const double tol = 1.e-14;
  const D_complex answer( -1.25, 0.25 );
  if ( abs( Z - answer ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
