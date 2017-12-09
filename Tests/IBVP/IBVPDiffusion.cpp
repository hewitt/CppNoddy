/// \file IBVPDiffusion.cpp
/// \ingroup Tests
/// \ingroup IBVP
/// Solving the heat diffusion equation
/// \f[ f_t = f_{yy} \f]
/// subject to \f$ f(0) = 0 \f$ and \f$ f(1) = 0 \f$
/// with initial condition \f$ f(y,t=0) = y(1-y) \f$.
/// The solution is computed over a range of \f$ t \f$
/// and the maximum deviation away from the analytical solution is found.
/// The test fails if this deviation is larger than a set tolerance \f$ \tau \f$.
/// The analytical solution is the Fourier series:
/// \f[ f(y,t) = \sum_{n=1,3,...}^M \frac{8}{n^3\pi^3} \sin ( n\pi y) \exp{ -n^2\pi^2 t } \f]
/// The series is truncated at the \f$ M \f$-th term, in such a way that the
/// first neglected term is of magnitude \f$ < \tau /10 \f$.

#include <IBVP_bundle.h>

// enumerate the system variables
enum { f, fd };

namespace CppNoddy
{
  namespace Example
  {
    class Diff_equation : public Equation_2matrix<double>
    {
    public:

      /// The problem is 2nd order and real
      Diff_equation() : Equation_2matrix<double> ( 2 ) {}

      /// Define the equation
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = 0.0;
      }

      /// Define the (BVP) deriv by providing the identity matrix
      void matrix0( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = 1.0;
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
      }
      
      /// To speed things up we'll overload this to say the mass matrix is constant
      void get_jacobian_of_matrix1_mult_vector( const DenseVector<double> &state, const DenseVector<double> &vec, DenseMatrix<double> &h  ) const
      {
        // blank definition leads to a zero result
      }
      
    };

    class Diff_both_BC : public Residual_with_coords<double>
    {
    public:
      // 1 BC and 2 unknowns + 1 time coordinate
      Diff_both_BC() : Residual_with_coords<double> ( 1, 2, 1 ) {}

      void residual_fn( const DenseVector<double>& z, DenseVector<double>& b ) const
      {
        b[ 0 ] = z[ f ];
      }
    };
  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== IBVP: An unsteady diffusion eqn =================\n";
  cout << "\n";

  // Diffusion equation
  Example::Diff_equation problem;
  // boundary conditions
  Example::Diff_both_BC BC_both;

  // domain definition
  double left = 0.0;
  double right = 1.0;
  // number of (spatial) nodal points
  unsigned ny = 201;
  // time step
  double dt = 0.005;
  // number of time steps
  unsigned max_steps = ( unsigned ) ( 3.0 / dt );
  // test tolerance
  double tol = 1.e-4;

  // construct our IBVP
  PDE_IBVP<double> heat( &problem, Utility::uniform_node_vector( left, right, ny ), &BC_both, &BC_both );

  for ( unsigned i = 0; i < ny; ++i )
  {
    double y = heat.solution().coord( i );
    heat.solution()( i, f ) = y * ( 1 - y );
    heat.solution()( i, fd ) = 1 - 2 * y;
  }

  // maximum difference between numerical and series solution at centre point
  double max_diff( 0.0 );
  // time step
  for ( unsigned i = 1; i < max_steps; ++i )
  {
    // take a time step
    try
    {
      heat.step2( dt );
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }

    // evaluate analytical Fourier series solution at y = 0.5
    unsigned mid = ( ny - 1 ) / 2;
    double y = left + mid * ( right - left ) / ( ny - 1 );
    double u( 0.0 );
    int en( -1 );
    double correction( 0.0 );
    do
    {
      en += 2;
      correction = 8 / ( std::pow( en * M_PI, 3 ) )
                   * std::exp( -std::pow( en * M_PI, 2 ) * heat.coord() ) * std::sin( en * M_PI * y );
      u += correction;
    }
    while ( std::abs( correction ) > tol / 10. );
    // examine the difference between the numerical and series solutions
    max_diff = std::max( max_diff, std::abs( u - heat.solution()( mid, f ) ) );
  }

  bool failed = true;
  if ( abs( max_diff ) < tol )
  {
    failed = false;
  }

  if ( failed )
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
