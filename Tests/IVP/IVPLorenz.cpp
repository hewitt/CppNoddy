/// \file IVPLorenz.cpp
/// \ingroup Tests
/// \ingroup IVP
/// Integrate the Lorenz equations
/// \f[ \dot x = a( y - x )\,,\quad \dot y = x ( b - z ) - y\,, \quad \dot z = -c z + xy\,, \f]
/// forward in time using an adaptive Runge-Kutta-Fehlberg routine. The
/// parameters are chosen to be \f$ a=10 \f$, \f$ b = 7 \f$, \f$ c = 8/3 \f$.
/// There is a fixed point solution \f$ (x,y,z) = (4,4,6) \f$, which is tested here.

#include <IVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Lorenz equations by inheriting Equation base class
    class Lorenz_equations : public Equation<double>
    {
    public:

      /// Construct a 3rd order ODE
      Lorenz_equations() : Equation<double>( 3 ) {};

      /// We implement the equation as 3 first-order ODEs
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[ 0 ] = a * ( z[ 1 ] - z[ 0 ] );
        f[ 1 ] = z[ 0 ] * ( b - z[ 2 ] ) - z[ 1 ];
        f[ 2 ] = -c * z[ 2 ] + z[ 0 ] * z[ 1 ];
      }

      /// The usual 3 parameters of the Lorenz eqns
      double a, b, c;

    };
  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== IVP: integrating the Lorenz system ==============\n";
  cout << "\n";
  cout << "   Running shoot method - fixed step choice \n";

  DenseVector<double> u( 3, 0.0 ), u_init( 3, 0.0 );
  // initialise
  u_init[ 0 ] = 1.0;
  u_init[ 1 ] = 1.0;
  u_init[ 2 ] = 1.0;
  const int num_of_steps( 5000000 );
  // set up the problem.
  Example::Lorenz_equations problem;
  problem.a = 10.0;
  problem.b = 7.0;
  problem.c = 8.0 / 3.0;

  // Construct an ODE from the problem.
  ODE_IVP<double> ode( &problem, 0.0, 200.0, num_of_steps );
  ode.store_every() = 10;

  const double tol = 1.e-7;
#ifdef TIME
  Timer timer;
  timer.start();
  do
  {
#endif
    try
    {
      u = ode.shoot45( u_init, tol, 0.1 );
    }
    catch (const std::runtime_error &error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }

#ifdef TIME
    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.stop();
  timer.print();
#endif

  cout.precision( 12 );
  if ( abs( u[ 0 ] - 4.0 ) > tol )
  {
    cout << "    Difference = " << abs( u[ 0 ] - 4.0 ) << "\n";
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
  ode.get_mesh().dump_gnu("./DATA/test.dat");
}
