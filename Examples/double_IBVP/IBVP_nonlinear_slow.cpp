/// \file IBVP_nonlinear_slow.cpp
/// \ingroup Examples
/// \ingroup IBVP_double
/// Solving the nonlinear problem
/// \f[ -u u_x - u_t + u_{yy} = y^2 t^2 e^{-2x} - y e^{-x}\,, \f]
/// subject to \f$ u(x=0,y,t) = ty \f$, \f$ u(x,t=0,y) = 0 \f$
/// and \f$ f(x, y = 0,t ) = 0 \f$, \f$ f(x, y = 1,t ) = t \f$.
/// The solution is computed over a range of \f$ t \f$ for \f$ (x,y) \in [0,10]\times [0,1] \f$
/// and the maximum deviation away from the exact solution \f$ u = yte^{-x} \f$ is found.
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
			return -y * std::exp( -x ) + std::pow( y * t * std::exp( -x ), 2 );
    }

    class nonlinear : public Equation_with_double_mass<double>
    {
    public:
      /// The problem is 2nd order and real
      nonlinear() : Equation_with_double_mass<double> ( 2 ) {}

      /// Define a nonlinear advection diffusion problem
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        // The system
        f[ U ] = z[ Ud ];
        f[ Ud ] = source( x(), y(), t() );
      }

      /// Define the unsteady terms by providing the mass matrix for x evolution
      void mass1( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // eqn 1 variable 0
        m( 1, 0 ) = -z[ U ];
      }

      /// Define the unsteady terms by providing the mass matrix for t evolution
      void mass2( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // eqn 1 variable 0
        m( 1, 0 ) = -1.0;
      }
            
    };

    // BOUNDARY CONDITIONS
    class BC_lower : public Residual_with_coords<double>
    {
    public:
    	// 1 constraint, 2nd order system, 2 coordinates (x & t)
      BC_lower() : Residual_with_coords<double> ( 1, 2, 2 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ U ];
      }
    };

    class BC_upper : public Residual_with_coords<double>
    {
    public:
    	// 1 constraint, 2nd order system, 2 coordinates (x & t)
      BC_upper() : Residual_with_coords<double> ( 1, 2, 2 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
      	const double x = coord( 0 );
      	const double t = coord( 1 );
        b[ 0 ] = z[ U ] - t * std::exp( - x );
      }
    };

  } // end Test namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== double_IBVP: non-constant mass matrix example ===\n";
  cout << "\n";

  // instantiate the problem
  Example::nonlinear problem;
  Example::BC_lower BC_bottom;
  Example::BC_upper BC_top;

  // domain definition
  double top = 1.0;
  double bottom = 0.0;
  double left = 0.0;
  double right = 10.0;
  // number of (spatial) nodal points
  unsigned ny = 11;
  unsigned nx = 1001;
  // time and time step
  double t_end = 1.;
  double dt = 0.01;

  // construct our IBVP
  PDE_double_IBVP<double> nlin( &problem, 
      Utility::uniform_node_vector( left, right, nx ), 
      Utility::uniform_node_vector( bottom, top, ny ), 
      &BC_bottom, &BC_top );

  // initial conditions 
  for ( unsigned i = 0; i < nx; ++i )
  {
    for ( unsigned j = 0; j < ny; ++j )
    {
      nlin.solution()( i, j, U ) = 0.0;
      nlin.solution()( i, j, Ud ) = 0.0;
    }
  }
  
  double max_error( 0.0 );
  int counter( 0 );
  do
  {
  	// inlet profile is time dependent, so we set it here for the next time level
    nlin.update_previous_solution();
    // now we can alter the solution stored in the double_IBVP class to define
    // the desired x=0 'initial' conditions at the next time step
    double t_next = nlin.t() + dt;
    for ( unsigned j = 0; j < ny; ++j )
    {
      double y = nlin.solution().coord( 0, j ).second;			
      nlin.solution()( 0, j, U ) = t_next * y;
      nlin.solution()( 0, j, Ud ) = t_next;			
    }
    nlin.step2( dt );
    for ( unsigned i = 0; i < nx; ++i )
	  {
      for ( unsigned j = 0; j < ny; ++j )
      {
        double x = nlin.solution().coord( i, j ).first;
        double y = nlin.solution().coord( i, j ).second;      
        double exact_u = y * nlin.t() * exp( - x );
        max_error = max( max_error, abs( exact_u - nlin.solution()( i, j, U ) ) );  
      }
	  }
    ++counter;
  } while ( nlin.t() < t_end ); 

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
