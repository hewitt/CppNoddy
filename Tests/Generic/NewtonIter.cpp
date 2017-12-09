/// \file NewtonIteration.cpp
/// \ingroup Tests
/// \ingroup Generic
/// A vector Newton iteration to find
/// a root of \f[ f(z) = z^3 - 1 \f] by splitting
/// it into a vector equation (real & imaginary parts).

#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Defines the problem: here z^3 -1 = 0, in vector form
    /// using real and imaginary parts.
    class VCube_root_problem : public Residual<double>
    {
    public:

      VCube_root_problem() : Residual<double>( 2 ) {}

      /// The residual function
      /// \param z The independent variable
      /// \param f The residual function f(z)
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[ 0 ] = z[ 0 ] * z[ 0 ] * z[ 0 ] - 3 * z[ 0 ] * z[ 1 ] * z[ 1 ] - 1.;
        f[ 1 ] = 3 * z[ 0 ] * z[ 0 ] * z[ 1 ] - z[ 1 ] * z[ 1 ] * z[ 1 ];
      }
    };
  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Newton: vector residual root test  ==============\n";
  cout << "\n";

  // Instantiate the problem
  Example::VCube_root_problem residual_problem;
  // A Newton object
  Newton<double> newton( &residual_problem );

  // Set an initial guess
  DenseVector<double> guess( 2, 0.0 );
  guess[ 0 ] = -0.5;
  guess[ 1 ] = 1.5;

  try
  {
    newton.iterate( guess );
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }

  const double tol = 1.e-7;
  if ( ( abs( guess[ 0 ] + 0.5 ) > tol ) ||
       ( abs( guess[ 1 ] - sqrt( 3. ) / 2. ) > tol ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    guess.dump();
    cout << abs( guess[ 1 ] - sqrt( 3. ) / 2. ) << "\n";
    cout << abs( guess[ 0 ] + 0.5 ) << "\n";
    return 1;
  }
  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;

}
