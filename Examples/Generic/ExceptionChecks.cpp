/// \file ExceptionChecks.cpp
/// \ingroup Examples
/// \ingroup Generic
/// We do some obviously dumb things to induce a handful
/// of common failures. The failures should result in exceptions
/// being thrown, which should be caught. The test is passed if
/// the required number of exceptions are successfully caught.

#include <Newton_bundle.h>
#include <DenseLinearSystem.h>

namespace CppNoddy
{
  namespace Example
  {
    /// The problem to be solved, here z^3 -1 in complex form.
    class Cube_root_problem : public Residual<D_complex>
    {
    public:

      Cube_root_problem() : Residual<D_complex>( 1 )
      {}

      /// The residual function for z^3 -1
      /// \param z The independent variable.
      /// \param f The value of the residual f(z).
      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &f ) const
      {
        f[ 0 ] = z[ 0 ] * z[ 0 ] * z[ 0 ] - 1.0;
      }
    };
  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Exception checks: testing forced failures =======\n";
  cout << "\n";


  int exceptions_caught( 0 );
  int tests( 0 );

  {
    ++tests;
    // set up noddy singular 2x2 problem
    DenseMatrix<double> A( 2, 2, 0.0 );
    DenseVector<double> B( 2, 0.0 );
    A( 0, 0 ) = 1.;
    A( 0, 1 ) = 2.;
    A( 1, 0 ) = 2.; // row 1 = 2 * row 0 -- so singular
    A( 1, 1 ) = 4.;
    B[ 0 ] = 5.;
    B[ 1 ] = 11.;

    DenseLinearSystem<double> system( &A, &B, "native" );

    // singular matrix should be caught in the native solver
    try
    {
      system.solve();
    }
    catch ( const std::exception & e )
    {
      ++exceptions_caught;
      cout << "  Caught a runtime exception in the native linear solver \n\n";
    }
  }

  {
    // too many potential geometry errors to test ... just the add method here
#ifdef PARANOID
    ++tests;
    // geometry checkes are only made under PARANOID flags
    DenseMatrix<double> C( 3, 3, 1.0 );
    DenseMatrix<double> D( 2, 2, 1.0 );
    try
    {
      C.add( D );
    }
    catch ( std::exception & e )
    {
      ++exceptions_caught;
      cout << "  Caught a geometry exception when adding two matrices. \n\n";
    }
#endif

  }

  {
    ++tests;
    // set the tolerance & max iteration number too low
    Example::Cube_root_problem problem;
    Newton<D_complex> newton( &problem, 8, 1.e-12 );
    D_complex half( 0.5, 0.5 );
    DenseVector<D_complex> guess( 1, half );
    try
    {
      newton.iterate( guess );
    }
    catch ( std::exception & e )
    {
      ++exceptions_caught;
      cout << "  Caught an iteration exception in a scalar Newton problem. \n\n";
    }
  }

  if ( tests != exceptions_caught )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << tests << " checks were run but only " << exceptions_caught
         << " exceptions were caught!\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
