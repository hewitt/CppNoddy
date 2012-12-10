/// \file IVP_Harmonic.cpp
/// \ingroup Examples
/// \ingroup IVP
/// Integrate the harmonic equation
/// \f[ f''(y) + \lambda f(y) = 0 \f] with \f$\lambda=10\f$
/// from \f$ y = 0\f$ to \f$1\f$ as an IVP,
/// using an adaptive Runge-Kutta-Fehlberg routine.

#include <cassert>

#include <Timer.h>
#include <IVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the harmonic equation by inheriting Equation base class
    class Harmonic_equation : public Equation<double>
    {
    public:
      /// The harmonic equation is 2nd order
      Harmonic_equation() : Equation<double> ( 2 )
      {
      }

      /// We implement the equation as 2 first-order ODEs
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[ 0 ] = z[ 1 ];
        f[ 1 ] = - lambda * z[ 0 ];
      }

      /// A parameter
      double lambda;

    };
  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{

  Timer timer;
  cout << "\n";
  cout << "=== IVP: solving harmonic equation ==================\n";
  cout << "   Solving between x = 0 to x = 10 \n";
  cout << "\n";
  cout << "   Running shoot45 method - auto step choice \n";

  // set up the problem
  Example::Harmonic_equation problem;
  problem.lambda = 10.0;

  // pass it to the ode
  ODE_IVP<double> ode( &problem, 0.0, 10.0, 1000 );
  // run the checks
  cout.precision( 12 );

  double tol( 1.e-7 );
  DenseVector<double> u_init( 2, 0.0 );
  DenseVector<double> u_final( 2, 0.0 );
  u_init[ 0 ] = 1.0;
  u_init[ 1 ] = 0.0;

  problem.coord(0) = 0.0;
  ode.shoot45( u_init, tol, 0.01 );

#ifdef TIME
  timer.start();
  do
  {
#endif
    try
    {
      u_final = ode.shoot45( u_init, tol, 0.01 );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    timer.counter()++;
#ifdef TIME
  }
  while ( timer.get_time() < 5000.0 );
  // stop and print the CPU time taken
  timer.stop();
  timer.print();
#endif

  bool failed( false );
  for ( int j = 1; j < 7; j++ )
  {
    // run through a few tolerance choices
    tol = 1. / ( pow( 10., j ) );
    // auto step choice; tolerance; initial step
    u_final = ode.shoot45( u_init, tol, 0.01 );
    // some output
    cout << "\n";
    cout << "    Error relative tol : " << tol << "\n";
    cout << "    Error |num - exact|: "
         << abs( u_final[ 0 ] - cos( sqrt( problem.lambda ) * 10. ) ) << "\n";

    if ( abs( u_final[ 0 ] - cos( sqrt( problem.lambda ) * 10. ) ) > 10. * tol )
    {
      failed = true;
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
