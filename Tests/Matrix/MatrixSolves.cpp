/// \file MatrixSolves.cpp
/// \ingroup Test
/// \ingroup Matrix
/// Example of the simple linear solvers using
/// a simple \f$ 2 \times 2 \f$ matrix problem.

#include <cassert>

#include <Timer.h>
#include <Types.h>
#include <Utility.h>
#include <DenseLinearSystem.h>
#include <BandedLinearSystem.h>

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Matrix: Example linear solver  ==================\n";
  cout << "\n";

  bool failed = false;
  // tolerance for the test
  const double tol = 1.e-10;

  //
  // SOLVE A SMALL "Dense"(!) 2X2 REAL SYSTEM with LAPACK or native
  //
  DenseMatrix<double> A( 2, 2, 0.0 );
  DenseVector<double> B( 2, 0.0 );
  A( 0, 0 ) = 1.;
  A( 0, 1 ) = 2.;
  A( 1, 0 ) = 3.;
  A( 1, 1 ) = 4.;
  B[ 0 ] = 5.;
  B[ 1 ] = 11.;

  std::cout << "     Simple 2x2 (dense) system solved by native Gauss Jordan routine.\n";
  DenseLinearSystem<double> small_system( &A, &B, "native" );

  try
  {
    small_system.solve();
  }
  catch ( const std::runtime_error &error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }
  DenseVector<double> answer( 2, 0.0 );
  answer[ 0 ] = 1.0;
  answer[ 1 ] = 2.0;
  B.sub( answer );
  if ( B.inf_norm() > tol )
  {
    std::cout << " Simple 2x2 system was not solved correctly\n";
    std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
    failed = true;
  }
  else
  {
    std::cout << "     Simple 2x2 `dense' solver works.\n";
  }

  //
  // CONCLUDING PASS/FAIL
  //
  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
}
