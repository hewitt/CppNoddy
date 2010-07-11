/// \file EVP_OrrSommerfeld_neutralcurve.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solves the following linear eigenvalue problem for values \f$ c \f$
/// that satisfy :
/// \f[ \phi''(y) - \alpha^2 \phi(y) - \psi(y) = 0\,, \f]
/// \f[ \psi''(y) - \alpha^2 \psi(y) - i \alpha Re \left \{ ( U(y) - c ) \psi(y) - U''(y) \phi \right \} = 0\,, \f]
/// subject to \f$ \phi(\pm 1) = \phi'(\pm 1) = 0 \f$ determining the
/// values of \f$ \alpha \f$ and \f$ Re \f$ that lead to \f$c_i = 0\f$.
/// In other words, we compute the temporal linear neutral curve
/// for perturbations to plane Poiseuille flow. The problem is solved by bolting
/// together the ODE_EVP (to get the initial global spectrum), the ODE_BVP
/// to (adaptively) refine the most dangerous eigenmode and the Newton class
/// to perform arc-length continuation of the neutral curve.

#include <cassert>

#include <EVP_bundle.h>
#include <Equation.h>
#include <Newton_bundle.h>


// enumerate the variables in the ODE
enum { phi, phid, psi, psid, eval };

namespace CppNoddy
{
  namespace Example
  {

    /// Globally define the Reynolds number and wavenumber
    double Re, alpha;
    /// The phase speed of the instability
    double wave_speed;
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
        // define the equation as four 1st order equations
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
      /// The OS local problem is a nonlinear 5th order complex BVP
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
        B[ 2 ] = z[ psi ] - 1.0; // an arbitrary amplitude traded against eigenvalue
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

    // Residual that defines the neutral curve
    class Neutral_residual : public Residual<double>
    {

    public:

      Neutral_residual( const unsigned& initial_N ) : Residual<double>( 1 )
      {
        // instantiate the object by doing a rough QZ pass of the eigenvalues
        OS_evp_equation evp;
        // instantiate the BCs
        OS_evp_both_BC BC_both;
        // domain is -1 to 1
        double left = -1.0;
        double right = 1.0;
        // number of nodal points
        int N = initial_N;
        // mesh
        DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
        // pass it to the ode
        ODE_EVP<D_complex> ode_global( &evp, nodes, &BC_both, &BC_both );
        try
        {
          // solve the global eigenvalue problem for growth rate
          ode_global.eigensolve();
        }
        catch ( std::runtime_error )
        {
          std::cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
          assert( false );
        }
        // get the eigenvalues near the origin
        ode_global.p_eigensystem() -> tag_eigenvalues_disc( + 1, 0.6 );
        // make a mesh of each tagged eigenmode in the ODE_EVP class
        ode_global.add_tagged_to_mesh();
        DenseVector<D_complex> lambdas = ode_global.p_eigensystem() -> get_tagged_eigenvalues();

        // instantiate the nonlinear BVP formulation
        bvp = new OS_bvp_equation;
        BC_left = new OS_bvp_left_BC;
        BC_right = new OS_bvp_right_BC;
        // pass it to the ode BVP solver
        ode_local = new ODE_BVP<D_complex>( bvp, nodes, BC_left, BC_right );
        ode_local -> set_monitor_det( false );
        // use our not-great initial guess from QZ code to start the solver off
        // here we just use the first (0) tagged eigenvalue
        ode_local -> solution() = ode_global.get_mesh( 0 );
        // interpolate onto a finer mesh
        ode_local -> solution().remesh1( Utility::uniform_node_vector( left, right, 201 ) );
        ode_local -> set_monitor_det( false );
        ode_local -> solve2();
        // adaptively refine the mesh
        for ( int i = 0; i < 8; ++i )
        {
          ode_local -> adapt( 1.e-2 );
          ode_local -> solve2();
          std::cout << "Adapted to " << ode_local -> solution().get_nnodes() << " nodes \n";
        }
      }

      ~Neutral_residual()
      {
        delete bvp;
        delete BC_left;
        delete BC_right;
        delete ode_local;
      }

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        Example::Re = z[ 0 ];
        // solve the nonlinear BVP
        ode_local -> solve2();
        // return the residual
        f[ 0 ] = ode_local -> solution()( 0, eval ).imag();
        wave_speed = ode_local -> solution()( 0, eval ).real();
      }

      ODE_BVP<D_complex>* ode_local;
      OS_bvp_equation* bvp;
      OS_bvp_left_BC* BC_left;
      OS_bvp_right_BC* BC_right;
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Temporal OSE via the equation interface ====\n";
  cout << "===  followed by local refinement of the eigenvalue  \n";
  cout << "===  and arclength continuation of the neutral curve.\n";
  cout << "\n";


  // screwed if we dont even have LAPACK
#ifndef LAPACK
  cout << " LAPACK support is a minumum requirement for this EVP computation.\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
  assert( false );
#endif

  double tol = 1.e-10;
  unsigned initial_QZ_mesh = 64;
  DenseVector<double> Re( 1, 5800.0 );
  Example::alpha = 0.99;
  Example::Re = Re[0];

  Example::Neutral_residual residual_problem( initial_QZ_mesh );

  Newton<double> newton( &residual_problem, 20, tol );
  newton.set_monitor_det( true );
  newton.rescale_theta() = true;
  newton.theta() = 0.0001;
  newton.init_arc( Re, &Example::alpha, 0.001, 0.01 );

  // set up the output file
  TrackerFile my_file( "DATA/neutral.dat" );
  my_file.precision( 8 );
  my_file.push_ptr( &Re, "Reynolds number" );
  my_file.push_ptr( &Example::alpha, "Wavenumber" );
  my_file.push_ptr( &Example::wave_speed, "Phase speed" );
  my_file.header();

  double min_Re( Re[ 0 ] );
  do
  {
    newton.arclength_solve( Re );
    my_file.update();
    min_Re = std::min( Re[ 0 ], min_Re );
  }
  while ( Re[ 0 ] < 6000.0 );

  std::cout << " Minimum Reynolds number (for this spatial resolution ) = " << min_Re << "\n";
  // the known minimum critical Re is 5772
  if ( std::abs( min_Re - 5772. ) < 0.5 )
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
  else
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }

}
