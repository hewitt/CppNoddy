/// \file IBVP_diffusion_nonlinear.cpp
/// \ingroup Examples
/// \ingroup IBVP
/// Solving the unstead diffusion problem:
/// \f[ U_t - U U_{yy} = -Ae^{-t} sin(\pi y) + ( (1+y) + Ae^{-t} sin(\pi y) )\pi^2 Ae^{-t} sin(\pi y)/\sigma \f]
/// with boundary conditions \f$ U(y=0)=1, U(y=1)=2 \f$; where \f$ A=10 \f$. The initial condition
/// is \f$ U(y,t=0)=1+y+A sin(\pi y) \f$.
/// The class constructs and solves the IBVP using 2nd order finite-differencing
/// in both space and time over a non-uniform spatial mesh. The output is compared to the exact
/// solution \f$ U(y,t) = =1+y+A sin(\pi y)e^{-t} \f$.

#include <cassert>

#include <PDE_IBVP.h>
#include <TrackerFile.h>
#include <OneD_Node_Mesh.h>
//
#include <Utility.h>
#include <Timer.h>

enum { U, Ud };

namespace CppNoddy
{
  namespace Example
  {

    double A( 10.0 );
    double sigma( 100.0 );

    class Diffusion_equations : public Equation_2matrix<double>
    {
    public:

      /// The problem is 2nd order and real
      Diffusion_equations() : Equation_2matrix<double> ( 2 ) {}

      /// Define the Karman equations
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      { 
        double E( A*exp( -coord(1) ) );
        double Y( coord(0) );
        // we will add some source terms to give an exact solution
        f[ 0 ] = z[ Ud ];
        f[ 1 ] = -E*sin(M_PI*Y) + ((1+Y)+E*sin(M_PI*Y))*(M_PI*M_PI)*E*sin(M_PI*Y)/sigma;
      }
      
      /// Define the domain-derivative terms by providing the mass matrix
      void matrix0( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // blank definition leads to a identity result
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = -z[ U ]/sigma;
      }

      /// To speed things up we'll overload this to say the mass matrix is constant
      void get_jacobian_of_matrix0_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

      /// Define the unsteady terms by providing the mass matrix
      void matrix1( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m( 1, 0 ) = 1.0;
      }

      /// To speed things up we'll overload this to say the mass matrix is constant
      void get_jacobian_of_matrix1_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }

    };

    /// Define the boundary conditions
    class Diffusion_left_BC : public Residual_with_coords<double>
    {
    public:
      // 3 BCs for 5 unknowns and one independent coordinate (time)
      Diffusion_left_BC() : Residual_with_coords<double> ( 1, 2, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ] - 1.0;
      }
    };

    class Diffusion_right_BC : public Residual_with_coords<double>
    {
    public:
      // 2 BCs for 5 unknowns and one independent coordinate (time)
      Diffusion_right_BC() : Residual_with_coords<double> ( 1, 2, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ] - 2.0;
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== IBVP: A nonlinear diffusion problem  ============\n";
  cout << "\n";

  // define the system
  Example::Diffusion_equations problem;
  // define the boundary conditions
  Example::Diffusion_left_BC BC_left;
  Example::Diffusion_right_BC BC_right;

  // domain definition
  double left = 0.0;
  double right = 1.0;
  // number of (spatial) nodal points
  unsigned ny = 400;
  // time step
  double dt = 0.005;
  // number of time steps
  unsigned max_steps = ( unsigned ) ( 2.0 / dt );

  // construct our IBVP with more nodes near the boundary
  PDE_IBVP<double> diffusion( &problem, Utility::power_node_vector( left, right, ny, 1.0 ), &BC_left, &BC_right );
  //
  for ( unsigned i = 0; i < ny; ++i )
  {
    double y( diffusion.solution().coord(i) );
    diffusion.solution()( i, U ) = (1+y) + Example::A*sin(M_PI*y);
    diffusion.solution()( i, Ud ) = 1 + Example::A*M_PI*cos(M_PI*y);
  }

  // set up output details
  TrackerFile metric( "./DATA/IBVP_diffusion_metric.dat" );
  metric.push_ptr( &diffusion.coord(), "time" );
  metric.push_ptr( &diffusion.solution()( 0, Ud ), "y=0 derivative" );
  metric.header();
  TrackerFile profs( "./DATA/IBVP_diffusion_profs.dat" );
  profs.push_ptr( &diffusion.coord(), "time" );
  profs.push_ptr( &diffusion.solution(), "solution" );
  profs.header();

  Timer timer;
  timer.start();
  //
  double max_error( 0.0 );
  // time step
  for ( unsigned i = 1; i < max_steps; ++i )
  {
    // take a time step
    try
    {
      diffusion.step2( dt );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    if ( i % 20 == 0 )
    {
      cout << diffusion.coord() << "\n";
      profs.update();
      profs.newline();
    }
    timer.counter()++;
    metric.update();
    //
    double yy( 0.5 );
    double Uexact = 1 + yy + Example::A*exp(-diffusion.coord())*sin(M_PI*yy);
    double error = diffusion.solution().get_interpolated_vars( 0.5 )[ U ] - Uexact;
    max_error = std::max( std::abs( error ), max_error );
    //
    // std::cout << " ****** " << diffusion.t() << " " << max_error << "\n";
  }
  timer.stop();
  timer.print();
  //
  const double tol( 1.e-4 );
  // check the BL transpiration vs the known solution
  if ( max_error > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " Max difference over time steps = " << max_error << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
