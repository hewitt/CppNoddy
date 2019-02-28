/// \file BVPKarmanAdaptive.cpp
/// \ingroup Tests
/// \ingroup BVP
/// Adaptively solve the Karman rotating-disk equations for the
/// flow above an infinite rotating disk:
/// \f[ U''(y) = U^2(y) + V(y)U'(y) - W^2(y) \f]
/// \f[ W''(y) = 2U(y)W(y) + V(y)W'(y)  \f]
/// \f[ 2U(y) + V'(y) = 0 \f]
/// with boundary conditions \f$ U(0)=V(0)=0 \f$, \f$ W(0)=1 \f$
/// and \f$ U(\infty ) \to 0 \f$, \f$ W(\infty ) \to 0 \f$.
/// The class constructs and solves the
/// global matrix problem using 2nd-order finite differences.

#include <BVP_bundle.h>

#include "../Utils_Fill.h"

// enumerate the 5 variables of the ODE system
enum {U, Ud, V, W, Wd};

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Karman equations
    class Karman_equations : public Equation_1matrix<double>
    {
    public:

      /// The Karman system is a 5th order real system of ODEs
      Karman_equations() : Equation_1matrix<double>( 5 ) {}

      /// Define the Karman system
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        // The 5th order system for ( U, U', V, W, W' )
        f[ U ] = z[ Ud ];
        f[ Ud ] = z[ U ] * z[ U ] + z[ V ] * z[ Ud ] - z[ W ] * z[ W ];
        f[ V ] = -2 * z[ U ];
        f[ W ] = z[ Wd ];
        f[ Wd ] = 2 * z[ U ] * z[ W ] + z[ V ] * z[ Wd ];
      }

      void matrix0( const DenseVector<double>&x, DenseMatrix<double> &m ) const
      {
        Utils_Fill::fill_identity(m);
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
        B[ 1 ] = z[ W ];
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: finite-difference soln of Karman eqns ======\n";
  cout << "===  Here we use an adaptive mesh after timing.   ===\n";
  cout << "\n";

  Example::Karman_equations problem;
  Example::Karman_left_BC BC_left;
  Example::Karman_right_BC BC_right;

  // Boundary layer is from 0 to 20
  double left = 0.0;
  double right = 40.0;
  // number of nodal points
  int N = 21;

  DenseVector<double> nodes = Utility::power_node_vector( left, right, N, 1.5 );

  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );

  for ( int i = 0; i < N; ++i )
  {
    double y = ode.solution().coord( i );
    ode.solution()( i, U ) = 0.0;
    ode.solution()( i, Ud ) = 0.0;
    ode.solution()( i, V ) = 0.0;
    ode.solution()( i, W ) = exp( -y );
    ode.solution()( i, Wd ) = -exp( -y );
  }

  try
  {
    ode.solve2();
  }
  catch (const std::runtime_error &error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }
  std::cout << "===  Adapting the mesh.\n";

  // output to check the adapted mesh
  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  TrackerFile my_file( "./DATA/BVP_Karman.dat", 12 );
  my_file.push_ptr( &ode.solution(), "soln" );

  const double tol( 1.e-4 );
  int adapt_counter ( 1 );
  do
  {
    // one adapt
    ode.adapt( 1.e-4 );
    ++adapt_counter;
    // solve
    ode.solve2();
    // no. of nodes in the adapted mesh
    N = ode.solution().get_nnodes();
    cout << " Adapted mesh to " << ode.solution().get_nnodes() << " nodes.\n";
    cout << " Adapted error = " << abs( ode.solution()( N - 1, V ) + 0.88447 ) << "\n";
    // adapt until condition is satisfied or max adaptions are done
  }
  while ( ( abs( ode.solution()( N - 1, V ) + 0.88447 ) > tol ) && ( adapt_counter < 10 ) );
  // store solution for plotting
  my_file.update();

  // check the BL transpiration vs the known solution
  if ( abs( ode.solution()( N - 1, V ) + 0.88447 ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " Difference = " << abs( ode.solution()( N - 1, V ) + 0.88447 ) << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

}
