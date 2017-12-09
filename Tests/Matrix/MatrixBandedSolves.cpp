/// \file MatrixBandedSolves.cpp
/// \ingroup Test
/// \ingroup Matrix
/// Example of a simple "banded" inear solver using
/// \f$ 2 \times 2 \f$ matrix problem is solved.
/// Then a penta-diagonal problem is solved using the
/// dense, banded and sparse containers. The native linear Gaussian
/// elimination solvers are used unless the PETSC_D/Z compiler
/// options are used, in which case the linear solver phase
/// calls the PETSc library as appropriate.

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
  cout << "=== Matrix: Example linear banded solver  ============\n";
  cout << "\n";

  bool failed = false;
  // tolerance for the test
  const double tol = 1.e-10;

  //
  // SOLVE A BANDED REAL SYSTEM using LAPACK or native
  //
  const unsigned offdiag = 2;
  // N = size of some of the larger tests
  const unsigned N = 511;
  const double D = 12 * ( 1. / ( N - 1 ) ) * ( 1. / ( N - 1 ) );
  DenseMatrix<double> AD( N, N, 0.0 );
  DenseVector<double> BD( N, D );
  BandedMatrix<double> AB( N, offdiag, 0.0 );
  DenseVector<double> BB( N, D );
  Utility::fill_band( AD, 0, -30.0 );
  Utility::fill_band( AD, -1, 16.0 );
  Utility::fill_band( AD, 1, 16.0 );
  Utility::fill_band( AD, -2, -1.0 );
  Utility::fill_band( AD, 2, -1.0 );
  Utility::fill_band( AB, 0, -30.0 );
  Utility::fill_band( AB, -1, 16.0 );
  Utility::fill_band( AB, 1, 16.0 );
  Utility::fill_band( AB, -2, -1.0 );
  Utility::fill_band( AB, 2, -1.0 );

  cout << "     Using dense matrix solver : " << N << "x" << N << " system \n";
  cout << "     Using the native dense routine\n";
  DenseLinearSystem<double> dense_system( &AD, &BD, "native" );

  try
  {
    dense_system.solve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }

  cout << "  * Not checking.\n";
  cout << "\n";
  cout << "     Comparing the banded matrix solver solution : ";
  cout << N << "x" << 2 * AB.noffdiag() + 1 << " system \n";

  cout << "     Using the native banded routine\n";
  BandedLinearSystem<double> banded_system( &AB, &BB, "native" );

  try
  {
    banded_system.solve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }

  BB.sub( BD );
  if ( std::abs( BB.two_norm() ) > tol )
  {
    failed = true;
    cout << " \033[1;31;48m * Banded solver does not give same result as dense solver \033[0m\n";
  }
  else
  {
    cout << "     Banded solver agrees with the dense solver.\n";
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
