/// \file BVP_Karman_arc.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Arc-length continuation of the Karman rotating-disk equations for the
/// flow above an infinite rotating disk:
/// \f[ U''(y) = U^2(y) + V(y)U'(y) - W^2(y) + s^2 \f]
/// \f[ W''(y) = 2U(y)W(y) + V(y)W'(y)  \f]
/// \f[ 2U(y) + V'(y) = 0 \f]
/// with boundary conditions \f$ U(0)=V(0)=0 \f$, \f$ W(0)=1 \f$
/// and \f$ U(\infty ) \to 0 \f$, \f$ W(\infty ) \to s \f$.
/// The class constructs and solves the
/// global matrix problem using 2nd-order finite differences and
/// employs arc-length continuation to obtain the first three solution
/// branches in the neighbourhood of \f$ s=0 \f$. The validation of the
/// results is based upon the approximate location of the second limit
/// point.

#include <cassert>

#include <BVP_bundle.h>
#include <Utility.h>
#include <TrackerFile.h>
#include <Timer.h>

// enumerate the 5 variables of the ODE system
enum {U, Ud, V, W, Wd};

namespace CppNoddy
{
  namespace Example
  {
    /// relative rotation rate
    double s;

    /// Define the Karman equations
    class Karman_equations : public Equation<double>
    {
    public:

      /// The Karman system is a 5th order real system of ODEs
      Karman_equations() : Equation<double>( 5 ) {}

      /// Define the Karman system
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        // The 5th order system for ( U, U', V, W, W' )
        f[ U ] = z[ Ud ];
        f[ Ud ] = z[ U ] * z[ U ] + z[ V ] * z[ Ud ] - z[ W ] * z[ W ] + s * s;
        f[ V ] = -2 * z[ U ];
        f[ W ] = z[ Wd ];
        f[ Wd ] = 2 * z[ U ] * z[ W ] + z[ V ] * z[ Wd ];
      }
    };

    /// Define the boundary conditions
    class Karman_left_BC : public Residual<double>
    {
    public:
      // 3 BCs for 5 unknowns
      Karman_left_BC() : Residual<double> ( 3, 5 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ U ];
        B[ 1 ] = z[ V ];
        B[ 2 ] = z[ W ] - 1.0;
      }
    };

    class Karman_right_BC : public Residual<double>
    {
    public:
      // 2 BCs for 5 unknowns
      Karman_right_BC() : Residual<double> ( 2, 5 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ U ];
        B[ 1 ] = z[ W ] - s;
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: Arc-length continuation of the Karman eqns ==\n";
  cout << "\n";

  Example::Karman_equations problem;
  Example::Karman_left_BC BC_left;
  Example::Karman_right_BC BC_right;

  // Boundary layer is from 0 to 20
  double left = 0.0;
  //need a large domain for the higher branch modes
  double right = 150.0;
  Example::s = 1.0;
  // number of nodal points
  unsigned N = 2001;
  // use a non-uniform mesh
  DenseVector<double> nodes = Utility::power_node_vector( left, right, N, 1.2 );

  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );

  // initial state is rigid-body rotation
  for ( unsigned i = 0; i < N; ++i )
  {
    ode.solution()( i, W ) = Example::s;
  }

  // output to check the adapted mesh
  TrackerFile my_file( "./DATA/BVP_Karman_arc.dat", 8 );
  my_file.push_ptr( &Example::s, "s" );
  my_file.push_ptr( &ode.solution()( N - 1, V ), "V_inf" );

  // initialise the arc-length routine
  double ds( -0.02 );
  ode.init_arc( &Example::s, ds, 0.01 );
  ode.rescale_theta() = true;
  ode.desired_arc_proportion() = 0.25;
  my_file.update();

  double last_approx_lp( 0.0 );
  // take 101 arc-length steps
  for ( int i = 0; i < 101; ++i )
  {
    try
    {
      ds = ode.arclength_solve( ds );
      //cout << Example::s << " " << ode.solution()( N-1, V ) << "\n";
    }
    catch ( ExceptionBifurcation )
    {
      cout << " Bifurcation detected near p = " << Example::s << "\n";
      last_approx_lp = Example::s;
    }
    my_file.update();
  }
  if ( std::abs( last_approx_lp + 0.1605 ) > 1.e-3 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " Difference = " << abs( last_approx_lp + 0.1605 ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
