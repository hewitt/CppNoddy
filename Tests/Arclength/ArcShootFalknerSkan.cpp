/// \file ArcShootFalknerSkan.cpp
/// \ingroup Tests
/// \ingroup Arclength
/// \ingroup BVP
/// Arc-length continue the Falkner-Skan equation
/// \f[ f'''(y) + f(y) f''(y) + \beta \left ( 1 - f'(y)^2 \right )= 0\,, \f]
/// for varying values of the
/// Hartree parameter \f$ \beta \f$ -- around the well known limit point. The boundary conditions are
/// \f$ f(0)=f'(0)=0 \f$ and \f$ f'(\infty) = 1 \f$. We start from the Blasius solution,
/// then arc-length continue to find the limit point at \f$ \beta \approx -0.19 \f$.
///
/// Using Runge-Kutta to solve a BVP is not a great method, but here it's used to
/// demo the arc-length continuation for scalar problems.

#include <cassert>

#include <IVP_bundle.h>
#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Falkner-Skan equation
    class FS_eqn : public Equation<double>
    {
    public:
      /// The FSE is a 3rd order ODE
      FS_eqn() : Equation<double>( 3 ) {}

      /// We implement the equation as 3 first-order real ODEs.
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        f[ 0 ] = z[ 1 ];
        f[ 1 ] = z[ 2 ];
        f[ 2 ] = -z[ 0 ] * z[ 2 ] - beta * ( 1.0 - z[ 1 ] * z[ 1 ] );
      }

      /// The Hartree parameter
      double beta;

    };

    /// Define a residual function using the boundary
    /// conditions for the Blasius profile.
    class FS_residual : public Residual<double>
    {
    public:
      ODE_IVP<double> *ode;
      FS_eqn *eqn;

      FS_residual() : Residual<double> ( 1 )
      {
        eqn = new FS_eqn;
        eqn -> beta = 0.0; // init the Hartree parameter
        ode = new ODE_IVP<double> ( eqn, 0.0, 20.0, 2000 );
      }

      ~FS_residual()
      {
        delete ode;
        delete eqn;
      }

      /// implement a residual function.
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        DenseVector<double> u( 3, 0.0 ); // 3 variables of FS problem
        u[ 0 ] = 0.0;          // f = 0
        u[ 1 ] = 0.0;          // f' = 0
        u[ 2 ] = z[ 0 ];       // f'' = unknown
        // shoot the ODE
        u = ode -> shoot4( u );
        // return the residual of f'(infty) = 1 condition
        f[ 0 ] = u[ 1 ] - 1.;
      }
    };
  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== ARC: arc-length continuation of the FS equation =\n";
  cout << "\n";

  // the residual function to be satisfied
  Example::FS_residual problem;

  // initial conditions
  problem.eqn -> beta = 0.0;
  DenseVector<double> stress( 1, 0.4 );
  // Scalar Newton iteration problem
  Newton<double> newton( &problem );

  newton.set_monitor_det( true );
  // initialise a state for arc length continuation
  newton.init_arc( stress, &problem.eqn -> beta, -0.01, 0.1 );

  // output for the data
  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  TrackerFile my_file( "./DATA/FSarc_path.dat" );
  my_file.push_ptr( &problem.eqn -> beta, "Hartree parameter" );
  my_file.push_ptr( &stress, "Shear stress at the wall" );
  my_file.header();

  // seek out the limit point in an ad hoc way
  double approx_limit_point = 0.0;
  do
  {
    double last_beta = problem.eqn -> beta;
    try
    {
      try
      {
        newton.arclength_solve( stress );
        //cout << problem.eqn -> beta << " " << stress[0] << "\n";
      }
      catch ( const std::runtime_error &error )
      {
        cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
        return 1;
      }
    }
    catch ( const ExceptionBifurcation &bifn )
    {
      cout << " Bifurcation detected between beta = " << last_beta
           << " and beta = " << problem.eqn -> beta << "\n";
      cout << " Continuing further.\n";
      newton.set_monitor_det( false );
      approx_limit_point = 0.5 * ( problem.eqn -> beta + last_beta );
    }
    my_file.update();
  }
  while ( std::abs( problem.eqn -> beta ) > 0.01 );

  if ( abs( approx_limit_point + 0.1988 ) > 0.0005 )
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
