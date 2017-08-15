/// \file Arc_circle_vector.cpp
/// \ingroup Examples
/// \ingroup Arclength
/// A simple arc-length continuation solving the vector equation
/// \f[ x^2 + p^2 = 1\,, y = \sin(x) \f] where \f$p\f$ is a parameter.
/// The solution is obtained by Newton iteration combined with arclength
/// continuation with the parameter \f$p\f$.

#include <cassert>

#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    class Arc_problem : public Residual<double>
    {
    public:
      // the control parameter
      double p;

      Arc_problem() : Residual<double>( 2 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[ 0 ] = 1.0 - std::pow( z[ 0 ], 2 ) - std::pow( p, 2 );
        f[ 1 ] = z[ 1 ] - std::sin( z[ 0 ] );
      }
    };
  }
}

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== ARC: cont of x^2 + p^2 = 1.0, y = sin(x) ========\n";
  cout << "\n";

  // initial guess
  DenseVector<double> x( 2, 0.0 );
  x[ 0 ] = 0.4;
  x[ 1 ] = 0.4;

  // Instantiate the problem
  Example::Arc_problem residual;
  // initialise the parameter
  residual.p = 0.9;
  // A Newton object
  Newton<double> newton( &residual );
  newton.set_monitor_det( true );
  newton.init_arc( x, &residual.p, 0.05, 0.1 );

  double tol( 1.e-7 );
  bool failed = false;
  for ( int i = 0; i < 400; ++i )
  {
    try
    {
      try
      {
        newton.arclength_solve( x );
      }
      catch ( std::runtime_error )
      {
        cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
        assert( false );
      }
    }
    catch ( ExceptionBifurcation )
    {
      cout << " Bifurcation detected near x = " << x[ 0 ]
           << " p = " << residual.p << "\n\n";
    }
    if ( abs( pow( x[ 0 ], 2 )
              + pow( residual.p, 2 ) - 1.0 ) > tol )
    {
      failed = true;
      cout << " Error = " << pow( x[ 0 ], 2 )
           + pow( residual.p, 2 ) - 1.0 << "\n";
    }
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}


