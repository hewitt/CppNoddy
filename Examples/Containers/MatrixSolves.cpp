/// \file MatrixSolves.cpp
/// \ingroup Examples
/// \ingroup Containers
/// Example of the simple linear solvers implemented
/// for dense, banded and sparse matrix objects. To begin
/// with a simple \f$ 2 \times 2 \f$ matrix problem is solved.
/// Then a penta-diagonal problem is solved using the
/// dense, banded and sparse containers. The native linear Gaussian
/// elimination solvers are used unless the LAPACK/SUPERLU compiler
/// options are used, in which case the linear solver phase
/// calls the LAPACK/SUPERLU library.

#include <cassert>

#include <Timer.h>
#include <Types.h>
#include <Utility.h>
#include <DenseLinearSystem.h>
#include <BandedLinearSystem.h>
#include <SparseLinearSystem.h>

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Matrix: Example linear solver  ==================\n";
  cout << "\n";

  //
  // SOLVE A SMALL 2X2 REAL SYSTEM
  //
  DenseMatrix<double> A( 2, 2, 0.0 );
  DenseVector<double> B( 2, 0.0 );
  A( 0, 0 ) = 1.;
  A( 0, 1 ) = 2.;
  A( 1, 0 ) = 3.;
  A( 1, 1 ) = 4.;
  B[ 0 ] = 5.;
  B[ 1 ] = 11.;

#ifdef LAPACK

  std::cout << " Simple 2x2 system solved by LAPACK LU.\n";
  DenseLinearSystem<double> small_system( &A, &B, "lapack" );
#else

  std::cout << " Simple 2x2 (dense) system solved by native Gauss Jordan routine.\n";
  DenseLinearSystem<double> small_system( &A, &B, "native" );
#endif

  try
  {
    small_system.solve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }
  DenseVector<double> answer( 2, 0.0 );
  answer[ 0 ] = 1.0;
  answer[ 1 ] = 2.0;
  B.sub( answer );
  const double tol = 1.e-10;
  bool failed = false;
  if ( B.inf_norm() > tol )
  {
    std::cout << " Simple 2x2 system was not solved correctly\n";
    std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
    failed = true;
  }
  else
  {
    std::cout << " Simple solver works.\n";
  }

  //
  // SOLVE A BANDED REAL SYSTEM
  //
  std::cout << "\n Penta-diagonal system check:";
  const unsigned N = 511;
  const unsigned offdiag = 2;
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
  
  Timer timer;
  
  cout << " Using dense matrix solver : " << N << "x" << N << " system \n";
  DenseMatrix<double> aD( AD );
  DenseVector<double> bD( BD );
#ifdef LAPACK

  cout << " Using the LAPACK LU dense routine\n";
  DenseLinearSystem<double> dense_system( &aD, &bD, "lapack" );
#else

  cout << " Using the native Gauss-Jordan dense routine\n";
  DenseLinearSystem<double> dense_system( &aD, &bD, "native" );
#endif

#ifdef TIME
  do
  {
#endif
    // reset the matrix
    aD = AD;
    bD = BD;
#ifdef TIME
    timer.start();
#endif
    try
    {
      dense_system.solve();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
#ifdef TIME
    timer.stop();
    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.print();
  timer.reset();
#endif

  cout << "  * Not checking.\n";
  cout << "\n";
  cout << " Comparing the banded matrix solver solution : ";
  cout << N << "x" << 2 * AB.noffdiag() + 1 << " system \n";

  BandedMatrix<double> aB( AB );
  DenseVector<double> bB( BB );
#ifdef LAPACK

  cout << " Using the LAPACK LU banded routine\n";
  BandedLinearSystem<double> banded_system( &aB, &bB, "lapack" );
#else

  cout << " Using the native Gauss-Jordan banded routine\n";
  BandedLinearSystem<double> banded_system( &aB, &bB, "native" );
#endif

#ifdef TIME
  do
  {
#endif
    // reset the matrix
    aB = AB;
    bB = BB;
#ifdef TIME
    timer.start();
#endif
    try
    {
      banded_system.solve();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
#ifdef TIME
    timer.stop();
    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.print();
  timer.reset();
#endif

  bB.sub( bD );
  if ( std::abs( bB.two_norm() ) > tol )
  {
    failed = true;
    cout << " \033[1;31;48m * Banded solver does not give same result as dense solver \033[0m\n";
  }

  cout << "\n";
  cout << " Comparing the sparse matrix solver solution : ";
  SparseMatrix<double> AS( N, N );
  Utility::fill_band( AS, 0, -30.0 );
  Utility::fill_band( AS, -1, 16.0 );
  Utility::fill_band( AS, 1, 16.0 );
  Utility::fill_band( AS, -2, -1.0 );
  Utility::fill_band( AS, 2, -1.0 );

  cout << N << " rows and " << AS.nelts() << " elts \n";
  cout << " Using the native Gauss-Jordan sparse routine\n";

  DenseVector<double> BS( N, D );
  SparseMatrix<double> aS( AS );
  DenseVector<double> bS( BS );

#ifdef LAPACK

  cout << " Using the SUPERLU sparse routine\n";
  SparseLinearSystem<double> sparse_system( &aS, &bS, "superlu" );
#else

  cout << " Using the native Gauss-Jordan sparse routine\n";
  SparseLinearSystem<double> sparse_system( &aS, &bS, "native" );
#endif

#ifdef TIME
  do
  {
#endif
    // reset the matrix
    aS = AS;
    bS = BS;
#ifdef TIME
    timer.start();
#endif
    try
    {
      sparse_system.solve();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
#ifdef TIME
    timer.stop();
    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.print();
  timer.reset();
#endif

  bS.sub( bD );

  if ( std::abs( bS.two_norm() ) > tol )
  {
    cout << " \033[1;31;48m * Sparse solver does not give same result as dense solver \033[0m\n";
    failed = true;
  }

  if ( failed )
  {
    cout << " || dense - sparse ||_2 = " << std::abs( bS.two_norm() ) << "\n";
    cout << " || dense - banded ||_2 = " << std::abs( bB.two_norm() ) << "\n";
  }

  //
  // CHECK COMPLEX SOLVERS FOR A SIMPLE SMALL SYSTEM
  //
  DenseMatrix<D_complex> CA( 2, 2, 0.0 );
  CA( 0, 0 ) = 1.;
  CA( 0, 1 ) = D_complex(2.,0.);
  CA( 1, 0 ) = 3.;
  CA( 1, 1 ) = 4.;
  DenseVector<D_complex> CB( 2, 0.0 );
  CB[ 0 ] = 5;
  CB[ 1 ] = 11;

#ifdef LAPACK

  std::cout << "\n Simple 2x2 system solved by complex LAPACK LU.\n";
  DenseLinearSystem<D_complex> small_Csystem( &CA, &CB, "lapack" );
#else

  std::cout << "\n Simple 2x2 (dense) complex system solved by native Gauss Jordan routine.\n";
  DenseLinearSystem<D_complex> small_Csystem( &CA, &CB, "native" );
#endif

  try
  {
    small_Csystem.solve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }
  CB.sub( answer );
  if ( CB.inf_norm() > tol )
  {
    std::cout << " Simple 2x2 COMPLEX system was not solved correctly\n";
    std::cout << " residual vector's inf_norm = " << CB.inf_norm() << "\n";
    failed = true;
  }
  else
  {
    std::cout << " Simple complex solver works.\n";
  }

  SparseMatrix<D_complex> CSA( 2, 2 );
  CSA( 0, 0 ) = 1.;
  CSA( 0, 1 ) = D_complex(2.,0.);
  CSA( 1, 0 ) = 3.;
  CSA( 1, 1 ) = 4.;
  CB[ 0 ] = 5;
  CB[ 1 ] = 11;

#ifdef SUPERLU

  std::cout << "\n Simple 2x2 system solved by complex SUPERLU.\n";
  SparseLinearSystem<D_complex> small_CSsystem( &CSA, &CB, "superlu" );
#else

  std::cout << "\n Simple 2x2 ('sparse') complex system solved by native Gauss Jordan routine.\n";
  SparseLinearSystem<D_complex> small_CSsystem( &CSA, &CB, "native" );
#endif
  try
  {
    small_CSsystem.solve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  CB.sub( answer );
  if ( CB.inf_norm() > tol )
  {
    std::cout << " Simple 2x2 COMPLEX 'sparse' system was not solved correctly\n";
    std::cout << " residual vector's inf_norm = " << CB.inf_norm() << "\n";
    failed = true;
  }
  else
  {
    std::cout << " Simple complex sparse solver works.\n";
  }

  //
  // CONCLUDING PASS/FAIL
  //
  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }


}
