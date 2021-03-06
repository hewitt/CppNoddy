/// \file EVPHarmonicEasy_lapack.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the harmonic equation
/// \f[ f''(x) + \lambda f(x) = 0 \f]
/// as an eigenvalue problem for \f$ \lambda \f$ over the unit domain
/// with homogeneous boundary conditions for \f$ f(x) \f$, returning
/// the smallest eigenvalue. This is the 'easy' approach, using the
/// ODE_EVP class which is constructed from an Equation_with_mass object, with
/// Residual objects for the boundary conditions. The solver then
/// only needs a pointer to the eigenvalue variable.

#include <EVP_bundle.h>
#include <Utility.h>
#include "../Utils_Fill.h"

// enumerate the variables in the ODE
enum {f, fd };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the harmonic equation by inheriting the Equation base class
    class harmonic_equation : public Equation_2matrix<double>
    {
    public:
      /// The harmonic equation is a 2nd order real ODE
      harmonic_equation() : Equation_2matrix<double>( 2 ) {}

      /// The harmonic equation
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = 0.0;
      }

      /// matrix to multiply the BVP coordinate
      void matrix0( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        Utils_Fill::fill_identity(m);
      }

      /// Define the eigenvalue terms by providing the mass matrix
      /// This defines the term lambda * z[ f ] ;
      void matrix1( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        // eigenvalue multiplies unknown 0 in equation 1
        m( 1, 0 ) = 1.0;
      }

    };

    class harmonic_both_BC : public Residual<double>
    {
    public:
      // 1 reisudal and 2 unknowns
      harmonic_both_BC() : Residual<double> ( 1, 2 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
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
  cout << "=== EVP: The harmonic equation done the easy way ====\n";
  cout << "\n";

  // test pass/fail
  bool failed( true );


  Example::harmonic_equation problem;
  Example::harmonic_both_BC BC_both;
  // domain is 0 to 1
  double left = 0.0;
  double right = 1.0;
  // number of nodal points
  int N = 256;
  // EV is pi^2 so we'll guess = 10 + 0i
  D_complex guess( 10.0, 0.0 );
  // pass/fail tolerance
  const double tol = 1.e-3;

  cout << " Solving the system using LAPACK:\n";
  // set up the ode eigenproblem
  ODE_EVP<double> ode_lapack( &problem, Utility::uniform_node_vector( left, right, N ), &BC_both, &BC_both );
  ode_lapack.p_eigensystem() -> set_calc_eigenvectors( true );
  ode_lapack.p_eigensystem() -> set_shift( guess );
  try
  {
    // solve the global eigenvalue
    ode_lapack.eigensolve();
  }
  catch (const std::runtime_error &error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }
  // get the eigenvalues in a disk of radius 1 about the guess
  ode_lapack.p_eigensystem() -> tag_eigenvalues_disc( + 1, 1. );
  // get the tagged eigenvalue(s) -- hopefully only 1
  DenseVector<D_complex> lapack_lambdas = ode_lapack.p_eigensystem() -> get_tagged_eigenvalues();
  if ( abs( lapack_lambdas[ 0 ].real() - M_PI * M_PI ) < tol )
  {
    failed = false;
    cout << "  LAPACK solver works.\n";
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }

  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;

}
