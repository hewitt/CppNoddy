/// \file EVP_trivially_linearised.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solving the model eigenvalue problem
/// \f[ - \sigma Re f(y) + f''(y) = - Re f(y) \f]
/// subject to \f$ f(0) = 0 \f$ and \f$ f(\pi) = 0 \f$
/// for the eigenvalue \f$\sigma \f$.
/// In this case we imagine that the above linearised
/// problem arises from considering the temporal evolution
/// of linear perturbations to steady solutions of the nonlinear problem
/// \f[ - F_t(y,t) + F_{yy}(y,t) = - Re \sin ( F(y,t) ) \f]
/// subject to \f$ F(0) = 0 \f$ and \f$ F(\pi) = 0 \f$.
/// The (in this case trivial) base state is computed using
/// the ODE_BVP class, but we also pass the mass matrix that
/// defines the linearised eigenproblem, thereby allowing us
/// to directly construct the generalised matrix problem for the
/// eigenvalues. The results are compared to the analytical
/// solution that \f$ \sigma = 1 - n^2/Re \f$ for $\f n = 1,2,3... \f$.

#include <cassert>

#include <BVP_bundle.h>
#include <EVP_bundle.h>

// enumerate the variables in the ODE
enum {f, fd };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Berman equation by inheriting the Equation base class
    class Model_equation : public Equation_with_mass<double>
    {
    public:

      /// A 2nd order real ODE
      Model_equation() : Equation_with_mass<double>( 2 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = -Re * std::sin( z[ f ] );
      }

      /// The mass matrix for the undstady term
      void mass( const DenseVector<double>& z, DenseMatrix<double>& m ) const
      {
        m( 1, f ) = -Re;
      }

      // The Reynolds number
      double Re;
    };

    class Model_BC : public Residual<double>
    {
    public:
      // one condition for a system with 2 unknowns      
      Model_BC() : Residual<double> ( 1, 2 ) {}

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
  cout << "=== EVP: a linearised EVP on top of a nonlinear BVP =\n";
  cout << "\n";

#ifndef LAPACK
  cout << " BLAS/LAPACK support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
#else

  Example::Model_equation problem;
  Example::Model_BC BC;

  // set the Reynolds number
  problem.Re = 0.5;
  // domain is -1 to 1
  double left = 0.0;
  double right = M_PI;
  // number of nodal points
  int N = 201;
  // mesh
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
  // error in validation check
  double max_error( 0. );
  // tolerance of validation check
  double tol( 1.e-4 );
  
  // pass the ode to the BVP class
  ODE_BVP<double> ode( &problem, nodes, &BC, &BC );
  ode.set_monitor_det( false );

  // a vector for the eigenvalues
  DenseVector<D_complex> sigmas;

  TrackerFile evals( "./DATA/EVP_evals.dat" );
  evals.push_ptr( &problem.Re, "Fake Reynolds number" );
  evals.push_ptr( &sigmas, "The complex eigenvalues" );
  do
  {  
    try
    {
      // solve for the (trivial) base flow.
      ode.solve2();
      // make 2 dense matrices
      DenseMatrix<double> a, b;
      // ask the ODE_BVP class to assemble the linearised EVP  
      ode.assemble_linear_evp( a, b );  
      // construct and solve the EVP
      DenseLinearEigenSystem<double> system( &a, &b );
      system.eigensolve();
      // find any eigenvalues in a disc around the origin
      system.tag_eigenvalues_disc( + 1, 10. );
      sigmas = system.get_tagged_eigenvalues();
      // dump the eigenvalues to disk for fun
      evals.update();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    // find the most unstable mode for validation purposes
    double most_unstable( sigmas[ 0 ].real() );
    for ( unsigned i = 1; i < sigmas.size(); ++i )
    {
      if ( sigmas[ i ].real() > most_unstable )
      {
        most_unstable = sigmas[ i ].real();
      }
    }
    max_error = std::max( std::abs( most_unstable - ( 1.0 - 1.0 / problem.Re ) ), max_error );
    problem.Re += 0.1;
  } while ( problem.Re < 3.0 );

  if ( std::abs( max_error ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "    Max error = " << max_error << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

#endif    
}
