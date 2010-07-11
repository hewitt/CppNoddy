/// \file EVP_local_Harmonic.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solves the harmonic equation
/// \f[ f''(x) + \lambda f(x) = 0 \f]
/// as a LOCAL eigenvalue problem for \f$ \lambda \f$ over the unit domain
/// with homogeneous boundary conditions for \f$ f(x) \f$ but using
/// the nonlinear BVP solver to refine a guess at an eigenvalue.
/// For example, we would have run an (expensive) global eigenvalue
/// search (EVP_Harmonic) then refined the eigenvalue of interest in
/// the way shown in this example.
///
/// Note: This approach (as it stands) is not as efficient as it could
/// be since we include N extra degrees of freedom (where N is the number
/// of mesh points) into the problem for only one eigenvalue. This retains
/// the banded structure, but it would perhaps be better to solve using
/// a sparse structure, or an appropriate multi-step algorithm for the
/// for the banded matrix.

#include <EVP_bundle.h>
#include <Equation.h>

// enumerate the 3 variables
enum { f, fd, lambda };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the harmonic equation by inheriting Equation base class
    class Harmonic_equation : public Equation<D_complex>
    {
    public:

      /// The eqn is a *3rd* order complex nonlinear ODE
      /// because the eigenvalue is an unknown
      Harmonic_equation() : Equation<D_complex> ( 3 ) {}

      /// Define the harmonic eqn
      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = - z[ lambda ] * z[ f ];
        // inclusion of eigenvalue everywhere maintains the bandedness,
        // albeit at an obvious cost.
        g[ lambda ] = 0.0;
      }
    };

    class Harmonic_left_BC : public Residual<D_complex>
    {
    public:
      Harmonic_left_BC() : Residual<D_complex> ( 2, 3 ) {}

      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
      {
        B[ 0 ] = z[ f ];
        B[ 1 ] = z[ fd ] - 1.0; // arbitrary amplitude
      }
    };

    class Harmonic_right_BC : public Residual<D_complex>
    {
    public:
      Harmonic_right_BC() : Residual<D_complex> ( 1, 3 ) {}

      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
      {
        B[ 0 ] = z[ f ];
      }
    };


  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Local refinement of an eigenvalue  =========\n";
  cout << "\n";

  cout << " Number of points : |local eigenvalue - pi^2| \n";

  // set up the problem
  Example::Harmonic_equation problem;
  // set up BCs
  Example::Harmonic_left_BC BC_left;
  Example::Harmonic_right_BC BC_right;
  // set up the domain from 0 to 1
  double left = 0.0;
  double right = 1.0;
  // set our guess for the eigenvalue -- this would come
  // from the global method in a less trivial problem.
  D_complex eigenvalue_guess = 9.0;
  bool failed( true );
  // loop through some meshes with increasing numbers of points
  for ( unsigned i = 6; i <= 10; ++i )
  {
    unsigned N( 0 );
    N = unsigned( std::pow( 2.0, double( i ) ) );
    // mesh
    DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
    // construct the ODE_BVP object
    ODE_BVP<D_complex> ode( &problem, nodes, &BC_left, &BC_right );
    // makes no sense to monitor the determinant here
    ode.set_monitor_det( false );
    // run through the mesh & set an initial guess
    for ( unsigned i = 0; i < N; ++i )
    {
      double x = ode.solution().coord( i );
      ode.solution()( i, f ) = x * ( 1 - x );
      ode.solution()( i, fd ) = 1.0 - 2 * x;
      ode.solution()( i, lambda ) = eigenvalue_guess;
    }
    // solve for the eigenvalue/eigenfunction
    ode.solve2();
    // an error measure
    double abs_error( std::abs( ode.solution()( 1, lambda ) - M_PI * M_PI ) );
    // output for fun
    cout << "  " << N << " : " << abs_error << "\n";
    if ( abs_error < 1.e-4 )
      failed = false;
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
