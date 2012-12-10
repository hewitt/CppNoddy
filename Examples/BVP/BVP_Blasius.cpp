/// \file BVP_Blasius.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Solving the Blasius equation
/// \f[ f'''(y) + f(y) f''(y) = 0\,, \f]
/// with \f$ f(0)=f'(0)=0 \f$ and \f$ f'(\infty) = 1 \f$
/// in the domain \f$ y \in [0,\infty ] \f$
/// by applying the ODE_BVP class.
/// The class constructs and solves the
/// global matrix problem using 2nd-order finite differences
/// over a finite domain of specified size.

#include <cassert>

#include <BVP_bundle.h>
#include <Equation_1matrix.h>
#include <Residual_with_coords.h>

// enumerate the 3 variables
enum { f, fd, fdd };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Blasius equation by inheriting Equation base class
    class Blasius_equation : public Equation_1matrix<double>
    {
    public:

      /// The Blasius eqn is a 3rd order real ODE
      Blasius_equation() : Equation_1matrix<double> ( 3 ) {}

      /// Define the Blasius eqn
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = z[ fdd ];
        g[ fdd ] = -z[ f ] * z[ fdd ];
      }
      
      void matrix0( const DenseVector<double>&x, DenseMatrix<double> &m ) const
      {
        Utility::fill_identity(m);
      }
    };

    class Blasius_left_BC : public Residual<double>
    {
    public:
      // 2 residuals and 3 unknowns
      Blasius_left_BC() : Residual<double> ( 2, 3 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ f ];
        B[ 1 ] = z[ fd ];
      }
    };

    class Blasius_right_BC : public Residual<double>
    {
    public:
      // 1 residual and 3 unknowns
      Blasius_right_BC() : Residual<double> ( 1, 3 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ fd ] - 1.0;
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: finite-difference solution of Blasius ======\n";
  cout << "\n";

  cout << " Number of points : approx. error \n";

  // equation
  Example::Blasius_equation problem;
  // boundary conditions
  Example::Blasius_left_BC BC_left;
  Example::Blasius_right_BC BC_right;

  double left = 0.0;        // from x=0
  double right = 20.0;      // to x=20
  bool failed = false;

  const double tol = 1.e-5; // pass/fail tolerance
  // let's repeat the computation for a sequence of N's
  // to check on convergence
  for ( int N = 32; N <= 2048; N *= 2 )
  {
    // mesh
    DenseVector<double> nodes( Utility::power_node_vector( left, right, N, 1.2 ) );
    // pass it to the ode
    ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );
    for ( int i = 0; i < N; ++i )
    {
      double y = ode.solution().coord( i );
      ode.solution()( i, f ) = y * ( 1.0 - exp( -y ) );
      ode.solution()( i, fd ) = ( 1.0 - exp( -y ) ) + y * exp( -y );
      ode.solution()( i, fdd ) = 2.0 * exp( -y ) - y * exp( -y );
    }

    try
    {
      ode.solve2();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }

    const double c = 1.65519036023e0;
    const double answer = 1. / pow( c, 1.5 );
    if ( abs( ode.solution()( 0, fdd ) - answer ) > tol )
    {
      failed = true;
    }
    else
    {
      failed = false;
    }

    // compare to the known stress value
    std::cout << "  " << N << " " << abs( ode.solution()( 0, fdd ) - answer ) << "\n";
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
