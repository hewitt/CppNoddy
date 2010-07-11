/// \file IBVP_linear.cpp
/// \ingroup Examples
/// \ingroup IBVP_double
/// Solving the linear equation
/// \f[ - u_x - u_t + u_{yy} = ( ( 1 - y^2 )( x + t ) - 2 ) e^{-xt}\,, \f]
/// subject to \f$ u(x=0,y,t) = 1-y^2 \f$ \f$ u(x,t=0,y) = 1 - y^2 \f$
/// and \f$ f(x, y = \pm 1,t ) = 0 \f$.
/// The solution is computed over a range of \f$ t \f$ for \f$ (x,y) \in [0,10]\times [-1,1] \f$
/// and the maximum deviation away from the exact solution \f$ u = (1-y^2) e^{-xt} \f$ is found.
/// The test fails if this deviation is larger than a set tolerance \f$ \tau \f$.
/// The example is solved using the PDE_double_IBVP for problems that are 
/// parabolic in 2 coordinates \f$ (x,t) \f$ with a BVP in \f$ y \f$.

#include <IBVP_double_bundle.h>
#include <Utility.h>

enum { U, Ud };

namespace CppNoddy
{
  namespace Example
  {

    // A source term rigged up to give a nice exact solution
    double source( const double& x, const double& y, const double& t )
    {
      return - 2 * std::exp( -x*t ) + ( 1 - y * y ) * ( x + t ) * std::exp( -x*t );
    }

    class diffusion_double : public Equation_with_double_mass<double>
    {
    public:
      /// The problem is 2nd order and real
      diffusion_double() : Equation_with_double_mass<double> ( 2 ) {}

      /// Define a nonlinear advection diffusion problem
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        // The system
        f[ U ] = z[ Ud ];
        f[ Ud ] = source( x(), y(), t() );
      }

      /// Provide the exact Jacobian rather than using finite-differences
      void jacobian( const DenseVector<double> &z, DenseMatrix<double> &jac ) const
      {
        jac( 0, Ud ) = 1.0;
      }

      /// Define the unsteady terms by providing the mass matrix
      void mass1( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // eqn 1 variable 0
        m( 1, 0 ) = -1.0;
      }

      void mass2( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // eqn 1 variable 0
        m( 1, 0 ) = -1.0;
      }
      
      void get_jacobian_of_mass1_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // constant mass1 matrix
      }
      void get_jacobian_of_mass2_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // constant mass2 matrix
      }

    };

    // BOUNDARY CONDITIONS
    class dd_BC : public Residual_with_coords<double>
    {
    public:
    	// 1 constraint, 2nd order system, 2 coordinates (x & t)
      dd_BC() : Residual_with_coords<double> ( 1, 2, 2 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ];
      }
    };

  } // end Test namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== double_IBVP: linear diffusive problem ===========\n";
  cout << "\n";

  // instantiate the problem
  Example::diffusion_double problem;
  Example::dd_BC BC_bottom;
  Example::dd_BC BC_top;

  // domain definition
  double top = 1.0;
  double bottom = -1.0;
  double left = 0.0;
  double right = 10.0;
  // number of (spatial) nodal points
  unsigned ny = 301;
  unsigned nx = 301;
  // time and time step
  double t_end = 1.0;
  double dt = 0.01;

  // construct our IBVP
  PDE_double_IBVP<double> diffusion( &problem, 
      Utility::uniform_node_vector( left, right, nx ), 
      Utility::uniform_node_vector( bottom, top, ny ), 
      &BC_bottom, &BC_top );
  
  // initial conditions 
  for ( unsigned i = 0; i < nx; ++i )
  {
    for ( unsigned j = 0; j < ny; ++j )
    {
      double y = diffusion.solution().coord( i, j ).second;
      diffusion.solution()( i, j, U ) = ( 1.0 - y * y );
      diffusion.solution()( i, j, Ud ) = - 2 * y ;
    }
  }
  
  double max_error( 0.0 );
  do
  {
    diffusion.update_previous_solution();
    diffusion.step2( dt );
    for ( unsigned i = 0; i < nx; ++i )
    {
      for ( unsigned j = 0; j < ny; ++j )
      {
        double x = diffusion.solution().coord( i, j ).first;
        double y = diffusion.solution().coord( i, j ).second;      
        double exact_u = ( 1 - y * y ) * exp( - x * diffusion.t() );
        max_error = max( max_error, abs( exact_u - diffusion.solution()( i, j, U ) ) );  
      }
    }
  } while ( diffusion.t() < t_end ); 

  const double tol( 1.e-4 );
  // check the BL transpiration vs the known solution
  if ( max_error > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " Difference = " << max_error << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
  
}
