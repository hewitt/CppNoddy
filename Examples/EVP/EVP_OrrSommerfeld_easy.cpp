/// \file EVP_OrrSommerfeld_easy.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solves the following linear eigenvalue problem for values \f$ c \f$
/// that satisfy :
/// \f[ \phi''(y) - \alpha^2 \phi(y) - \psi(y) = 0\,, \f]
/// \f[ \psi''(y) - \alpha^2 \psi(y) - i \alpha Re \left \{ ( U(y) - c ) \psi(y) - U''(y) \phi \right \} = 0\,, \f]
/// subject to \f$ \phi(\pm 1) = \phi'(\pm 1) = 0 \f$ where
/// \f$ \alpha = 1.02 \f$, \f$ Re = 5772.2 \f$ and \f$ U(y) = 1 - y^2 \f$.
/// These values approximately correspond to the first neutral temporal mode
/// in plane Poiseuille flow, therefore the test to be satisfied is that an eigenvalue
/// exists with \f$ \Im ( c ) \approx 0 \f$. Here we use the ODE_EVP class with
/// a very low number of mesh points to get an approximation to the temporal
/// spectrum of modes. The 'most dangerous' mode is then extracted from this
/// set, interpolated onto a MUCH finer mesh, and used as the initial guess
/// in a nonlinear boundary-value problem using the ODE_BVP class. This leads
/// to a local refinement of the selected mode that could equally be applied
/// to other eigenvalues very simply.

#include <cassert>

#include <Types.h>
#include <Equation.h>
#include <EVP_bundle.h>

// enumerate the variables in the ODE
enum { phi, phid, psi, psid, eval };

namespace CppNoddy
{
  namespace Example
  {

    /// Globally define the Reynolds number and wavenumber
    double Re, alpha;
    /// Globally define the base flow
    double U( double y )
    {
      return 1.0 - y * y;
    };
    /// Globally define the base flow curvature
    double Udd( double y )
    {
      return -2.0;
    };

    /// Define the OS equation for the global QZ EVP
    class OS_evp_equation : public Equation_with_mass<D_complex>
    {
    public:
      /// The OS equation is a 4th order complex ODE
      OS_evp_equation() : Equation_with_mass<D_complex>( 4 ) {}

      /// The OS equation
      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &g ) const
      {
        // define the equation as 4 1st order equations
        g[ phi ] = z[ phid ];
        g[ phid ] = z[ psi ] + alpha * alpha * z[ phi ];
        g[ psi ] = z[ psid ];
        g[ psid ] = alpha * alpha * z[ psi ]
                    + D_complex( 0.0, 1.0 ) * alpha * Re * ( U( y() ) * z[ psi ] - Udd( y() ) * z[ phi ] );
      }

      /// Define the unsteady terms by providing the mass matrix
      /// Here we define the eigenvalue contribution to the g[ psid ] equation
      ///   - D_complex( 0.0, 1.0 ) * alpha * Re * c * z[ psi ]
      void mass( const DenseVector<D_complex>& z, DenseMatrix<D_complex>& m ) const
      {
        // the eigenvalue is in equation 3, and multiplies unknown 2
        m( 3, 2 ) = D_complex( 0.0, 1.0 ) * alpha * Re;
      }
    };


    /// Define the OSE for the local refinement procedure
    class OS_bvp_equation : public Equation<D_complex>
    {
    public:
      /// The OS *LOCAL* problem is a nonlinear 5th order complex BVP
      /// because the eigenvalue is added to the 4th order equation
      OS_bvp_equation() : Equation<D_complex>( 5 ) {}

      /// The OS equation
      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &g ) const
      {
        // define the equation as 5 1st order equations
        g[ phi ] = z[ phid ];
        g[ phid ] = z[ psi ] + alpha * alpha * z[ phi ];
        g[ psi ] = z[ psid ];
        g[ psid ] = alpha * alpha * z[ psi ]
                    + D_complex( 0.0, 1.0 ) * alpha * Re * ( U( y() ) * z[ psi ] - Udd( y() ) * z[ phi ] )
                    - D_complex( 0.0, 1.0 ) * alpha * Re * z[ eval ] * z[ psi ];
        g[ eval ] = 0.0;
      }

    };

    // BOUNDARY CONDITIONS FOR THE EVP SOLVER - same at both boundaries
    class OS_evp_both_BC : public Residual<D_complex>
    {
    public:
      // 2 boundary conditions and 4 unknowns
      OS_evp_both_BC() : Residual<D_complex> ( 2, 4 ) {}

      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
      {
        B[ 0 ] = z[ phi ];
        B[ 1 ] = z[ phid ];
      }
    };


    // BOUNDARY CONDITIONS FOR THE BVP SOLVER
    class OS_bvp_left_BC : public Residual<D_complex>
    {
    public:
      // 3 boundary conditions and 5 unknowns
      OS_bvp_left_BC() : Residual<D_complex> ( 3, 5 ) {}

      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
      {
        B[ 0 ] = z[ phi ];
        B[ 1 ] = z[ phid ];
        B[ 2 ] = z[ psi ] - 1.0; // an arbitrary amplitude traded against the eigenvalue
      }
    };

    class OS_bvp_right_BC : public Residual<D_complex>
    {
    public:
      // 2 boundary conditions and 5 unknowns
      OS_bvp_right_BC() : Residual<D_complex> ( 2, 5 ) {}

      void residual_fn( const DenseVector<D_complex> &z, DenseVector<D_complex> &B ) const
      {
        B[ 0 ] = z[ phi ];
        B[ 1 ] = z[ phid ];
      }
    };


  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Temporal OSE via the equation interface ====\n";
  cout << "===  followed by local refinement of the eigenvalue.\n";
  cout << "\n";


  // screwed if we dont even have LAPACK
#ifndef LAPACK
  cout << " LAPACK support is a minumum requirement for this EVP computation.\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
  assert( false );
#endif

  // instantiate the eigenproblem
  Example::OS_evp_equation evp;
  // instantiate the boundary conditions -- same applies at both boundaries
  Example::OS_evp_both_BC BC_both;
  Example::Re = 5772.2;
  Example::alpha = 1.02;
  // domain is -1 to 1
  double left = -1.0;
  double right = 1.0;
  // number of nodal points
  int N = 32;
  // mesh
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );

  // pass it to the ode
  ODE_EVP<D_complex> ode_global( &evp, nodes, &BC_both, &BC_both );
  try
  {
    // solve the global eigenvalue problem on ropey mesh
    ode_global.eigensolve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }
  // get the eigenvalues near the origin
  ode_global.p_eigensystem() -> tag_eigenvalues_disc( + 1, 0.6 );
  // make a mesh of each tagged eigenmode in the ODE_EVP class
  ode_global.add_tagged_to_mesh();

  // instantiate the nonlinear BVP formulation
  Example::OS_bvp_equation bvp;
  // instantiate the BCs
  Example::OS_bvp_left_BC BC_left;
  Example::OS_bvp_right_BC BC_right;
  // pass it to the ode BVP solver
  ODE_BVP<D_complex> ode_local( &bvp, nodes, &BC_left, &BC_right );
  // use our not-great initial guess from QZ code to start the solver off
  // here we just use the first (0) tagged eigenvalue
  ode_local.solution() = ode_global.get_mesh( 0 );
  // interpolate onto a finer mesh
  ode_local.solution().remesh1( Utility::uniform_node_vector( left, right, 1024 ) );
  ode_local.set_monitor_det( false );
  // solve the nonlinear BVP
  ode_local.solve2();
  // check the eigenvalue
  const double tol = 1.e-4;
  bool failed( true );
  if ( abs( ode_local.solution()( 0, eval ).imag() ) < tol )
    failed = false;
  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
