/// \file Shoot_Berman.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Solving the fourth-order Berman (porous channel) similarity equation
/// \f[ f'''(y) = Re \left ( f(y)f''(y) - f'(y)^2 - K \right )\f]
/// where \f$y\in [-1,1]\f$ and \f$K\f$ is a pressure constant and \f$Re\f$ is the
/// Reynolds number based on the channel half-height and the wall suction. The boundary
/// conditions are that \f$f(\pm 1) = \pm 1\f$ and \f$f'(\pm 1 ) = 0\f$.
/// This BVP is solved via Runge-Kutta and (vector) Newton iteration. This similarity
/// form is an exact reduction of the Navier-Stokes equations.
///
/// Solving a BVP by converting it to an IVP like this is generally not a smart move!

#include <cassert>

#include <Timer.h>
#include <IVP_bundle.h>
#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Berman equation by
    /// inheriting Equation base class.
    class Berman_equation : public Equation<double>
    {
    public:

      /// Construct the Berman equation as a 4th order real system
      Berman_equation() : Equation<double> ( 4 ) {}

      /// We implement the equation as 4 first-order ODEs.
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        f[ 0 ] = z[ 1 ];
        f[ 1 ] = z[ 2 ];
        f[ 2 ] = Re * ( - z[ 3 ] - z[ 1 ] * z[ 1 ] + z[ 0 ] * z[ 2 ] );
        f[ 3 ] = 0.0;
      }

      /// The Reynolds number
      double Re;

    };

    /// Define a residual function using the boundary
    /// conditions for the Berman equation.
    class Berman_residual : public Residual<double>
    {
    public:
      ODE_IVP<double>* ode;

      // There are two residuals to be satisfied
      Berman_residual() : Residual<double> ( 2 )
      {
        // instantiate the Berman equation as our problem
        eqn = new Berman_equation;
        // set the Reynolds number
        eqn -> Re = 1.0;
        // instantiate an ODE in [-1.,1] with a maximum of
        // 100 integrator steps, using this problem equation
        ode = new ODE_IVP<double>( eqn, -1.0, 1.0, 500 );

        // let's keep the history of the IVP for output in main
        ode -> set_store_soln( true );
      }

      ~Berman_residual()
      {
        delete ode;
        delete eqn;
      }

      /// implement a residual function.
      void residual_fn( const DenseVector<double>& unknown, DenseVector<double>& BC ) const
      {
        DenseVector<double> u( 4, 0.0 ); // 4 variables of Berman problem

        u[ 0 ] = -1.0;         // V = 0
        u[ 1 ] = 0.0;          // V' = 0
        u[ 2 ] = unknown[ 0 ]; // V'' = unknown0
        u[ 3 ] = unknown[ 1 ]; // A = unknown1

        u = ode -> shoot4( u );

        BC[ 0 ] = u[ 0 ] - 1.;
        BC[ 1 ] = u[ 1 ];
      }
    private:
      Berman_equation* eqn;
    };
  } //end Example namespace
} //end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== BVP_Shoot: vector Newton solution of Berman eqn =\n";
  cout << "\n";

  Example::Berman_residual problem;
  // instantiate a Vector Newton class for the Biharmonic
  Newton<double> Berman( &problem );

  DenseVector<double> guess( 2, 0.0 );
  guess[ 0 ] = -2.0;
  guess[ 1 ] = 1.0;

#ifdef TIME

  Timer timer;
  timer.start();
  do
  {
#endif
    try
    {
      Berman.iterate( guess );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }

#ifdef TIME
    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.stop();
  timer.print();
#endif

  // output the data
  TrackerFile my_file( "./DATA/Berman.dat" );
  my_file.push_ptr( &( problem.ode -> get_mesh() ) );
  my_file.update();

  // check against the known pressure constant K for Re=1.
  if ( std::abs( guess[ 1 ] - 0.70457 ) > 1.e-5 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout.precision( 14 );
    cout << guess[ 0 ] << " " << guess[ 1 ] << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }


}
