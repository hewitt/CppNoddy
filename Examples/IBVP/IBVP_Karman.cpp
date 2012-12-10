/// \file IBVP_Karman.cpp
/// \ingroup Examples
/// \ingroup IBVP
/// Solving the unsteady Karman rotating-disk equations for the
/// flow above an infinite rotating disk in a rotating fluid:
/// \f[ -U_t + U_{yy} = U^2 + VU_y - W^2 - W_\infty^2 \f]
/// \f[ -W_t + W_{yy} = 2U W + V W_y  \f]
/// \f[ 2U + V_y = 0 \f]
/// with boundary conditions \f$ U(y=0)=V(y=0)=0 \f$, \f$ W(y=0)=\Omega(t) \f$
/// and \f$ U(y \to \infty ) \to 0 \f$, \f$ W(y \to \infty ) \to 0 \f$. The
/// initial condition corresponds to rigid rotation of the fluid and boundary at
/// the frequency \f$ W_{\infty} \f$
/// The class constructs and solves the IBVP using 2nd order finite-differencing
/// in both space and time over a non-uniform spatial mesh.
/// The example shows the case \f$ W_\infty = 0 \f$ and
/// \f[ \Omega(t) = 1 - (1 - W_\infty) \exp( -10t )\f] and simply
/// checks for convergence to the correct steady state solution.

#include <cassert>

#include <Equation_2matrix.h>
#include <PDE_IBVP.h>
#include <TrackerFile.h>
#include <OneD_Node_Mesh.h>
//
#include <Utility.h>
#include <Timer.h>

enum { U, Ud, V, W, Wd };

namespace CppNoddy
{
  namespace Example
  {
    // rotation frequency at infinity, appears as a pressure gradient term
    double W_inf( 0.0 );
    // rotation of the disk, appears in the BCs
    double W_disk( 1.0 );

    class Karman_equations : public Equation_2matrix<double>
    {
    public:

      /// The problem is 5th order and real
      Karman_equations() : Equation_2matrix<double> ( 5 ) {}

      /// Define the Karman equations
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        // The 5th order system for ( U, U', V, W, W' )
        f[ 0 ] = z[ 1 ];
        f[ 1 ] = z[ 0 ] * z[ 0 ] + z[ 2 ] * z[ 1 ] - z[ 3 ] * z[ 3 ] + W_inf * W_inf;
        f[ 2 ] = -2 * z[ 0 ];
        f[ 3 ] = z[ 4 ];
        f[ 4 ] = 2 * z[ 0 ] * z[ 3 ] + z[ 2 ] * z[ 4 ];
      }
      
      /// Define the BVP deriv by providing the matrix
      void matrix0( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // identity matrix
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = 1.0;
        m( 2, 2 ) = 1.0;
        m( 3, 3 ) = 1.0;        
        m( 4, 4 ) = 1.0;        
      }

      /// To speed things up we'll overload this to say the mass matrix is constant
      void get_jacobian_of_matrix0_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

      /// Define the unsteady terms by providing the mass matrix
      void matrix1( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m( 1, 0 ) = -1.0;
        m( 4, 3 ) = -1.0;
      }
      
      /// To speed things up we'll overload this to say the mass matrix is constant
      void get_jacobian_of_matrix1_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

    };

    /// Define the boundary conditions
    class Karman_left_BC : public Residual_with_coords<double>
    {
    public:
      // 3 BCs for 5 unknowns and one independent coordinate (time)
      Karman_left_BC() : Residual_with_coords<double> ( 3, 5, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        // get the time level for unsteady boundary condition
        const double t = coord( 0 );
        const double W_bc = Example::W_disk + ( Example::W_inf - Example::W_disk ) * std::exp( - 2 * t * t );
        b[ 0 ] = z[ U ];
        b[ 1 ] = z[ V ];
        b[ 2 ] = z[ W ] - W_bc;
      }
    };

    class Karman_right_BC : public Residual_with_coords<double>
    {
    public:
      // 2 BCs for 5 unknowns and one independent coordinate (time)
      Karman_right_BC() : Residual_with_coords<double> ( 2, 5, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ];
        b[ 1 ] = z[ W ] - W_inf;
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== IBVP: The unsteady Karman equations  ============\n";
  cout << "\n";

  // define the system
  Example::Karman_equations problem;
  // define the boundary conditions
  Example::Karman_left_BC BC_left;
  Example::Karman_right_BC BC_right;

  // domain definition
  double left = 0.0;
  double right = 40.0;
  // number of (spatial) nodal points
  unsigned ny = 400;
  // time step
  double dt = 0.005;
  // number of time steps
  unsigned max_steps = ( unsigned ) ( 30.0 / dt );

  // construct our IBVP with more nodes near the boundary
  PDE_IBVP<double> karman( &problem, Utility::power_node_vector( left, right, ny, 2.0 ), &BC_left, &BC_right );
  //
  for ( unsigned i = 0; i < ny; ++i )
  {
    karman.solution()( i, U ) = 0.0;
    karman.solution()( i, Ud ) = 0.0;
    karman.solution()( i, V ) = 0.0;
    karman.solution()( i, W ) = Example::W_inf;
    karman.solution()( i, Wd ) = 0.0;
  }

  // set up output details
  TrackerFile metric( "./DATA/IBVP_Karman_metric.dat" );
  metric.push_ptr( &karman.coord(), "time" );
  metric.push_ptr( &karman.solution()( ny - 1, V ), "Ekman suction" );
  metric.header();
  TrackerFile profs( "./DATA/IBVP_Karman_profs.dat" );
  profs.push_ptr( &karman.coord(), "time" );
  profs.push_ptr( &karman.solution(), "solution" );
  profs.header();

  Timer timer;
  timer.start();
  // time step
  for ( unsigned i = 1; i < max_steps; ++i )
  {
    // take a time step
    try
    {
      karman.step2( dt );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    timer.counter()++;
    metric.update();
  }
  timer.stop();
  timer.print();

  const double tol( 1.e-4 );
  // check the BL transpiration vs the known solution
  if ( abs( karman.solution()( ny - 1, V ) + 0.88447 ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " Difference = " << abs( karman.solution()( ny - 1, V ) + 0.88447 ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
