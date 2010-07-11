/// \file Poisson_C.cpp
/// \ingroup Examples
/// \ingroup Poisson
/// Solving a Cartesian Poisson problem:
/// \f[ \nabla^2 \psi(x,y) = 2( x^2 + y^2 ) \f]
/// with \f[ \psi(x,\pm 1) = x^2\,, \quad\mbox{and}\quad \psi(\pm 1, y) = y^2 \f]
/// where \f$ (x,y) \in [-1,1]\times[-1,1] \f$.
/// The global problem is solved (one step) and result is compared to the
/// exact solution \f[ \psi(x,y) = x^2y^2 \f].

#include <cassert>

#include <Types.h>
#include <Poisson_Cartesian.h>
#include <Utility.h>
#include <Timer.h>

namespace CppNoddy
{
  namespace Example
  {

    /// Define the source function for the Poisson solver
    double source_fn( double &x, double &y )
    {
      return 2 * ( x * x + y * y ) ;
    }

    /// The function that defines the BCs
    double BC_fn( double &x, double &y )
    {
      return x * x * y * y;
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
  double tol = 1.e-9;
  // Dense matrices
  DenseMatrix<double> source( N, N, 1.0 );
  DenseMatrix<double> error( N, N, 1.0 );

  DenseVector<double> X = Utility::uniform_node_vector( -2.0, 2.0, N );
  DenseVector<double> Y = Utility::uniform_node_vector( -1.0, 1.0, N );

  Timer timer;
  cout << "\n";
  cout << "=== Poisson: Cartesian geometry =====================\n";
  cout << "    Solving in [-2 , 2] x [-1 , 1]      \n";
  cout << "    Using a " << N << "x" << N << " mesh\n";
  cout << "\n";

  // instantiate a Poisson_Cartesian problem
  Poisson_Cartesian problem( X[ 0 ], X[ N - 1 ], Y[ 0 ], Y[ N - 1 ], N, N, &source );

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
      // set up the boundary conditions
      for ( int i = 0; i < N; ++i )
      {
        for ( int j = 0; j < N; ++j )
        {
          if ( ( i == 0 ) || ( i == N - 1 ) || ( j == 0 ) || ( j == N - 1 ) )
          {
            // set boundary conditions
            source.set( i, j ) = Example::BC_fn( X[ i ], Y[ j ] );
          }
          else
          {
            source.set( i, j ) = Example::source_fn( X[ i ], Y[ j ] );
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
  } while ( timer.get_time() <  5000.0 );
  timer.stop();
  timer.print();
#endif

  // check the streamfunction -> velocity utility method too whilst we're here
  DenseMatrix<double> u( N, N, 0.0 );
  DenseMatrix<double> v( N, N, 0.0 );
  DenseMatrix<double> error_uv( N, N, 0.0 );
  Utility::vels_from_streamfn_Cartesian( source, X[1] - X[0], Y[1] - Y[0], u, v );
  for ( int i = 0; i < N; ++i )
  {
    for ( int j = 0; j < N; ++j )
    {
      // compute the errors
      error.set( i, j ) = source.get( i, j ) - pow( X[ i ], 2 ) * pow( Y[ j ], 2 );
      error_uv.set( i, j ) = std::max( u( i, j ) - 2 * X[ i ] * X[ i ] * Y[ j ],
            v( i, j ) + 2 * X[ i ] * Y[ j ] * Y[ j ] );
    }
  }


#ifdef DEBUG
  cout << "inf_norm of final error matrix  = " << error.inf_norm() << "\n";
  cout << "one_norm of final error matrix  = " << error.one_norm() << "\n";
  cout << "two_norm of final error matrix  = " << error.two_norm() << "\n";
  cout << "frob_norm of final error matrix = " << error.frob_norm() << "\n";
  cout << "inf_norm of velocity matrix = " << error_uv.inf_norm() << "\n";
#endif


  if ( std::max( error.inf_norm(), error_uv.inf_norm() ) > tol )
  {
    failed = true;
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "inf_norm of final error matrix  = " << error.inf_norm() << "\n";
    cout << "inf_norm of final error_uv matrix  = " << error_uv.inf_norm() << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
