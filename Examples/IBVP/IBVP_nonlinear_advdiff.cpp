/// \file IBVP_nonlinear_advdiff.cpp
/// \ingroup Examples
/// \ingroup IBVP
/// Solving the nonlinear advection diffusion equation
/// \f[ U_{yy} - Re U U_t = S( y, t ) \f]
/// subject to \f$ U(y=0,t) = 1 \f$ and \f$ U(y=1,t) = 2 \f$
/// with initial condition
/// \f$ U(y,t=0) = 1 + y + \epsilon y(1-y) \f$ and
/// a source term
/// \f[ S(y,t) = ( 1 + y + \epsilon y(1-y)e^{-t} ) \epsilon y(1-y ) e^{-t} - 2\epsilon e^{-t} / Re \f]
/// for some constant parameters \f$\epsilon \f$ and \f$ Re \f$.
/// This source term is chosen to allow an exact solution
/// \f[ U(y,t) = 1 + y + \epsilon y(1-y) e^{-t} \f]
/// against which the numerical computation is compared, by computing
/// the maximum absolute deviation over all spatial and temporal nodes.

#include <cassert>

#include <IBVP_bundle.h>

enum { U, Ud };

namespace CppNoddy
{
  namespace Example
  {

    // Global parameters
    // A diffusivity measure
    double Re;
    // Amplitude of a kick to the initial profile
    const double eps( 2.0 );
    // A source term rigged up to give a nice exact solution
    double source( const double& y, const double& t )
    {
      return ( 1 + y + eps * y * ( 1 - y ) * std::exp( -t ) ) * eps * y * ( 1 - y ) * std::exp( -t )
             - 2 * eps * std::exp( - t ) / Re;
    }

    class Nlin_adv_equation : public Equation_2matrix<double>
    {
    public:

      /// The problem is 2nd order and real
      Nlin_adv_equation() : Equation_2matrix<double> ( 2 ) {}

      /// Define a nonlinear advection diffusion problem
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        // The system
        f[ U ] = z[ Ud ];
        f[ Ud ] = source( coord(0), coord(1) );
      }

      /// Define the derivative terms by providing the mass matrix -- identity in this case
      void matrix0( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m(0,0)=1;m(1,1)=1;
      }

      void get_jacobian_of_matrix0_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        /// constant mass matrix, so we'll overlload this as empty to speed things up
      }

      /// Define the unsteady terms by providing the mass matrix
      void matrix1( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m( 1, 0 ) = - Re * z[ U ];
      }

      void get_jacobian_of_matrix1_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        h( 1, U ) = - Re * vec[ 0 ];
      }
      
    };

    // BOUNDARY CONDITIONS
    class Nlin_adv_left_BC : public Residual_with_coords<double>
    {
    public:
      Nlin_adv_left_BC() : Residual_with_coords<double> ( 1, 2, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ] - 1.0;
      }
    };

    class Nlin_adv_right_BC : public Residual_with_coords<double>
    {
    public:
      Nlin_adv_right_BC() : Residual_with_coords<double> ( 1, 2, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ] - 2.0;
      }
    };


  } // end Test namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== IBVP: Nonlinear advection diffusion =============\n";
  cout << "\n";

  // instantiate the problem
  Example::Nlin_adv_equation problem;
  Example::Nlin_adv_left_BC BC_left;
  Example::Nlin_adv_right_BC BC_right;

  // domain definition
  double left = 0.0;
  double right = 1.0;
  // number of (spatial) nodal points
  unsigned ny = 400;
  // number of time steps
  unsigned max_steps = 4000;
  // time step
  double dt = 1.0 / ( max_steps - 1 );

  // initial rotation frequency
  Example::Re = 1.0;

  // construct our IBVP
  PDE_IBVP<double> advection( &problem, Utility::uniform_node_vector( left, right, ny ), &BC_left, &BC_right );

  // initial conditions
  for ( unsigned i = 0; i < ny; ++i )
  {
    double y = advection.solution().coord( i );
    advection.solution()( i, U ) = 1.0 + y + Example::eps * y * ( 1 - y );
    advection.solution()( i, Ud ) = 1.0;
  }

  // set up output details
  TrackerFile my_file( "./DATA/IBVP_NlinAdvection.dat" );
  my_file.push_ptr( &advection.coord(), "time" );
  my_file.push_ptr( &advection.solution(), "Solution profile" );
  my_file.header();

  double max_err( 0.0 );
  // time step
  for ( unsigned i = 1; i < max_steps; ++i )
  {
    // take a time step
    try
    {
      advection.step2( dt );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    if ( i % 100 == 0 )
    {
      my_file.update();
      my_file.newline();
    }
    // generate the exact solution to compare to
    for ( unsigned i = 0; i < ny; ++i )
    {
      double y = advection.solution().coord( i );
      const double exact_solution = 1.0 + y + Example::eps * y * ( 1 - y ) * std::exp( -advection.coord() );
      max_err = std::max( max_err, std::abs( advection.solution()( i, 0 ) - exact_solution ) );
    }
  }

  if ( max_err > 1.e-6 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
#ifdef DEBUG
    cout.precision( 14 );
    cout << 1. / ( ny - 1 ) << " " << dt << " " << max_err << "\n";
#endif
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
