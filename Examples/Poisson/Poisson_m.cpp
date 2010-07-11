/// \file Poisson_m.cpp
/// \ingroup Examples
/// \ingroup Poisson
/// Solving a Poisson problem in the meridional plane of a cylinder:
/// \f[ \nabla^2 \psi(r,z) = 2 r^2 \f]
/// with \f[ \psi(r,\pm 1) = r^2\,, \quad\mbox{and}\quad \psi(0, z) = \psi(1,z) = 0 \f]
/// where \f$ (r,z) \in [0,1]\times[-1,1] \f$.
/// The global problem is solved (in one step) and result is compared to the
/// exact solution \f[ \psi(r,z) = r^2z^2 \f].

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
      return 4 * z * z + 2 * r * r;
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

  DenseVector<double> r = Utility::uniform_node_vector( 0.0, 2.0, N );
  DenseVector<double> z = Utility::uniform_node_vector( -1.0, 1.0, N );

  Timer timer;
  cout << "\n";
  cout << "=== Poisson: Meridional plane of a cylinder =========\n";
  cout << "    Solving in [0 , 2] x [-1 , 1] \n";
  cout << "    Using a " << N << "x" << N << " mesh\n";
  cout << "\n";

  // instantiate a Poisson_meridional problem
  Poisson_meridional problem( r[ 0 ], r[ N - 1 ], z[ 0 ], z[ N - 1 ], N, N, &source );

  cout << "    Global solver  : ";
#ifdef LAPACK

  cout << "  Using LAPACK banded solver.\n";
#else
  cout << "  Using native banded solver.\n";
#endif

#ifdef TIME
  timer.start();
  do
  {
#endif
    try
    {

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

      problem.solve();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    timer.counter()++;
#ifdef TIME
  }
  while ( timer.get_time() <  5000.0 );
  timer.stop();
  timer.print();
#endif

  for ( int i = 0; i < N; ++i )
  {
    for ( int j = 0; j < N; ++j )
    {
      // compute the error
      error.set( i, j ) = source.get( i, j ) - pow( r[ i ] * z[ j ], 2 );
    }
  }

#ifdef DEBUG
  cout << "inf_norm of final error matrix  = " << error.inf_norm() << "\n";
  cout << "one_norm of final error matrix  = " << error.one_norm() << "\n";
  cout << "two_norm of final error matrix  = " << error.two_norm() << "\n";
  cout << "frob_norm of final error matrix = " << error.frob_norm() << "\n";
#endif

  if ( error.inf_norm() > tol )
  {
    failed = true;
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "inf_norm of final error matrix  = " << error.inf_norm() << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
