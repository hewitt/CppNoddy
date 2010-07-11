/// \file Poisson_Stokes.cpp
/// \ingroup Examples
/// \ingroup Poisson
/// Solving a Poisson-like problem in the meridional plane of a cylinder
/// for the Stokes streamfunction
/// \f[ D^2 \psi(r,z) = 2 r^2 \f]
/// with \f[ \psi(r,\pm 1) = r^2\,, \quad\mbox{and}\quad \psi(0, z) = \psi(1,z) = 0 \f]
/// where \f$ (r,z) \in [0,1]\times[-1,1] \f$ and
/// \f[ D^2 \equiv \frac{\partial^2}{\partial r^2} - \frac{1}{r}\frac{\partial}{\partial r} +\frac{\partial^2}{\partial z^2} \f]
/// The global problem is solved (in one step) and result is compared to the
/// exact solution \f[ \psi(r,z) = r^2z^2 \f]

#include <cassert>

#include <Types.h>
#include <Poisson_meridional.h>
#include <Utility.h>
#include <Timer.h>

namespace CppNoddy
{
  namespace Example
  {

    /// Define the source function for the Poisson solver
    double source_fn( double &r, double &z )
    {
      return 2 * r * r;
    }

    /// The function that defines the BCs
    double BC_fn( double &r, double &z )
    {
      return r * r * z * z;
    }

  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  // Number of points in a square mesh.
  int N = 64;
  bool failed = false;
  double tol = 1.e-10;
  // Dense matrices
  DenseMatrix<double> source( N, N, 1.0 );
  DenseMatrix<double> error( N, N, 1.0 );

  DenseVector<double> r = Utility::uniform_node_vector( 0.0, 1.0, N );
  DenseVector<double> z = Utility::uniform_node_vector( 0.0, 1.0, N );

  Timer timer;
  cout << "\n";
  cout << "=== Poisson: Stokes streamfunction  =================\n";
  cout << "    Solving in [0 , 1] x [-1 , 1] \n";
  cout << "    Using a " << N << "x" << N << " mesh\n";
  cout << "\n";

  // instantiate a Poisson_meridional problem
  Poisson_meridional problem( r[ 0 ], r[ N - 1 ], z[ 0 ], z[ N - 1 ], N, N, &source );
  problem.set_stokes_streamfn();

  // set up the boundary conditions
  for ( int i = 0; i < N; ++i )
  {
    for ( int j = 0; j < N; ++j )
    {
      if ( ( i == 0 ) || ( i == N - 1 ) || ( j == 0 ) || ( j == N - 1 ) )
      {
        // set boundary conditions
        source.set( i, j ) = Example::BC_fn( r[ i ], z[ j ] );
      }
      else
      {
        source.set( i, j ) = Example::source_fn( r[ i ], z[ j ] );
      }
    }
  }


  cout << "    Global solver  : ";
#ifdef LAPACK
  cout << "  Using LAPACK banded solver.\n";
#else
  cout << "  Using native banded solver.\n";
#endif

  problem.solve();

  // check the streamfunction -> velocity utility method too whilst we're here
  DenseMatrix<double> u( N, N, 0.0 );
  DenseMatrix<double> w( N, N, 0.0 );
  DenseMatrix<double> error_uw( N, N, 0.0 );
  Utility::vels_from_streamfn_Stokes( source, r[1] - r[0], z[1] - z[0], u, w );
  for ( int i = 0; i < N; ++i )
  {
    for ( int j = 0; j < N; ++j )
    {
      // compute the errors
      error.set( i, j ) = source.get( i, j ) - pow( r[ i ] * z[ j ], 2 );
      error_uw.set( i, j ) = std::max( u( i, j ) + 2 * r[ i ] * z[ j ],
            w( i, j ) - 2 * z[ j ] * z[ j ] );
    }
  }

  // the velocity error probably will not converge with N ... because of the
  // 1/r behaviour & the streamfunction coming from a 2nd-order scheme.
  if ( ( error.inf_norm() > tol ) || ( error_uw.inf_norm() > 1.e-7 )  )
  {
    failed = true;
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "inf_norm of error matrix  = " << error.inf_norm() << "\n";
    cout << "inf_norm of error_uw matrix  = " << error_uw.inf_norm() << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
