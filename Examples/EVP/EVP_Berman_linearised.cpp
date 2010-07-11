/// \file EVP_Berman_linearised.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solving the model eigenvalue problem
/// \f[ - \sigma Re f''(y) + f^{(iv)}(y) = Re ( F(y) f'''(y) + f(y) F'''(y) - F'(y) f''(y) - f'(y) F''(y) ) \f]
/// subject to \f$ f(\pm 1 ) = f'(\pm 1 ) = 0 \f$ for the eigenvalue \f$\sigma \f$.
/// In this case we imagine that the above linearised
/// problem arises from considering the temporal evolution
/// of linear perturbations to steady solutions of the nonlinear problem
/// \f[ -F_t( y,t) + F^(iv)(y) = Re ( F(y)F'''(y) - F'(y)F''(y) ) \f]
/// subject to \f$ F(\pm 1 ) =  \pm 1 \f$ and \f$ F'(\pm 1 ) = 0 \f$
/// The (in this case trivial) base state is computed using
/// the ODE_BVP class, but we also pass the mass matrix that
/// defines the linearised eigenproblem, thereby allowing us
/// to directly construct the generalised matrix problem for the
/// eigenvalues. The validation check is that a neutral eigenvalue
/// exists near \f$ Re = 6\f$.

#include <cassert>
#include <BVP_bundle.h>
#include <EVP_bundle.h>
#include <TrackerFile.h>

// enumerate the variables in the ODE
enum {f, fd, fdd, fddd };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Berman equation by inheriting the Equation base class
    class Berman_equation : public Equation_with_mass<double>
    {
    public:

      /// The Berman equation is a 4th order real ODE
      Berman_equation() : Equation_with_mass<double>( 4 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = z[ fdd ];
        g[ fdd ] = z[ fddd ];
        g[ fddd ] = Re * ( z[ f ] * z[ fddd ] - z[ fd ] * z[ fdd ] );
      }

      /// The mass matrix for the undstady term
      void mass( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m( 3, fdd ) = -Re;
      }
      
      // The Reynolds number
      double Re;
    };

    class Berman_left_BC : public Residual<double>
    {
    public:
      // 2 boundary conditions and 4 unknowns
      Berman_left_BC() : Residual<double> ( 2, 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ f ] + 1.0;
        B[ 1 ] = z[ fd ];
      }
    };

    class Berman_right_BC : public Residual<double>
    {
    public:
      // 2 boundary conditions and 4 unknowns
      Berman_right_BC() : Residual<double> ( 2, 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ f ] - 1.0;
        B[ 1 ] = z[ fd ];
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: linear temporal stability of Berman eqn ====\n";
  cout << "\n";

#ifndef LAPACK
  cout << " BLAS/LAPACK support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
#else

  Example::Berman_equation problem;
  Example::Berman_left_BC BC_left;
  Example::Berman_right_BC BC_right;

  // set the Reynolds number
  problem.Re = 6.0;
  // domain is -1 to 1
  double left = -1.0;
  double right = 1.0;
  // number of nodal points
  int N = 301;
  // Example tolerance
  const double tol = 1.e-3;

  // mesh
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
  // pass it to the ode
  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );
  ode.set_monitor_det( false );
  
  // our initial guess
  for ( int i = 0; i < N; ++i )
  {
    double y = ode.solution().coord( i );
    // set f(y)
    ode.solution()( i, f ) = 1.5 * ( y - y * y * y / 3 );
    // set f'(y)
    ode.solution()( i, fd ) = 1.5 * ( 1 - y * y );
    // set f''(y)
    ode.solution()( i, fdd ) = -3 * y;
    // set f''(y)
    ode.solution()( i, fddd ) = -3;
  }
  
  DenseVector<D_complex > sigmas;
  // solve the problem using 2nd order finite-difference
  try
  {
    ode.solve2();
    // make 2 dense matrices
    DenseMatrix<double> a, b;
    // ask the ODE_BVP class to assemble the linearised EVP  
    ode.assemble_linear_evp( a, b );
    // construct and solve the EVP
    DenseLinearEigenSystem<double> system( &a, &b );
    system.eigensolve();
    // find any eigenvalues in a disc around the origin
    system.tag_eigenvalues_disc( + 1, 1. );
    sigmas = system.get_tagged_eigenvalues();  
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }
 
  if ( std::abs( sigmas.inf_norm() ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << std::abs( sigmas.inf_norm() ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

#endif

}
