/// \file ArcCircle.cpp
/// \ingroup Tests
/// \ingroup Arclength
/// A simple arc-length continuation solving the equation
/// \f[ x^2 + p^2 = 2\,, \f] where \f$p\f$ is a parameter. The solution is obtained by
/// Newton iteration combined with (pseudo) arclength continuation with
/// the parameter \f$p\f$. For a (slightly) more serious example see the
/// FalknerSkan_arc.cpp example.

#include <cassert>

#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the residual for arc-length continuation of a circle
    class Arc_problem : public Residual<double>
    {
    public:
      double p;

      Arc_problem() : Residual<double>( 1 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[ 0 ] = z[ 0 ] * z[ 0 ] + p * p - 2.0;
      }
    };

  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== ARC: continuation of x^2 + p^2 = 2.0 ============\n";
  cout << "\n";

  // Instantiate the problem
  Example::Arc_problem residual_problem;
  residual_problem.p = 1.0;

  const double tol = 1.e-10;
  // Scalar Newton iteration problem
  Newton<double> newton( &residual_problem, 8, tol );
  newton.set_monitor_det( true );

  // initial guess
  DenseVector<double> state( 1, 1.1 );
  // initialise a state for arc length cont.
  newton.init_arc( state, &residual_problem.p, 0.001, 0.1 );

  bool failed = false;
  for ( int i = 0; i < 100; ++i )
  {
    try
    {
      try
      {
        newton.arclength_solve( state );
      }
      catch ( std::runtime_error )
      {
        cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
        return 1;
      }
    }
    catch ( ExceptionBifurcation )
    {
      cout << " Bifurcation detected near x = " << state[ 0 ] << " p = " << residual_problem.p << "\n\n";
    }
    if ( abs( pow( state[ 0 ], 2 ) + pow( residual_problem.p, 2 ) - 2.0 ) > tol )
    {
      failed = true;
      cout << " Error = " << pow( state[ 0 ], 2 ) + pow( residual_problem.p, 2 ) - 2.0 << "\n";
    }
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
}
