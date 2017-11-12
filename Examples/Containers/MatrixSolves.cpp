/// \file MatrixSolves.cpp
/// \ingroup Examples
/// \ingroup Containers
/// Example of the simple linear solvers implemented
/// for dense, banded and sparse matrix objects. To begin
/// with a simple \f$ 2 \times 2 \f$ matrix problem is solved.
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
#include <SparseLinearSystem.h>

// #include "mpi.h"
// #include <PETSc.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  #if defined(PETSC_Z) || defined(PETSC_D)
  PetscInitialize(NULL,NULL,(char*)0,(char*)0);
  #endif

  cout << "\n";
  cout << "=== Matrix: Example linear solver  ==================\n";
  cout << "\n";

  bool failed = false;
  // tolerance for the test
  const double tol = 1.e-10;

  cout << "\n";
  cout << "   ================ SIMPLE 2X2 TEST =================\n";
  {
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
    if ( B.inf_norm() > tol )
    {
      std::cout << " Simple 2x2 system was not solved correctly\n";
      std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
      failed = true;
    }
    else
    {
      std::cout << " Simple 2x2 `dense' solver works.\n";
    }
  }


  cout << "\n";
  cout << "   ================ BANDED (PENTA TEST) =============\n";
  {
    //
    // SOLVE A BANDED REAL SYSTEM using LAPACK or native
    //
    std::cout << "\n Penta-diagonal system check:";
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
    timer.start();
    do
    {
    #endif
      // reset the matrix
      aD = AD;
      bD = BD;
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
      timer.counter()++;
    }
    while ( timer.get_time() < 5000.0 );
    timer.stop();
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
    timer.start();
    do
    {
    #endif
      // reset the matrix
      aB = AB;
      bB = BB;
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
      timer.counter()++;
    }
    while ( timer.get_time() < 5000.0 );
    timer.stop();
    timer.print();
    timer.reset();
    #endif

    bB.sub( bD );
    if ( std::abs( bB.two_norm() ) > tol )
    {
      failed = true;
      cout << " \033[1;31;48m * Banded solver does not give same result as dense solver \033[0m\n";
    }
    else
    {
      cout << " Banded solver agrees with the dense solver.\n";
    }

    {
      //
      // SOLVE the BANDED REAL SYSTEM as above but as a sparse system using native and superlu
      //
      cout << "\n";
      cout << " Comparing the sparse matrix solver solution : ";
      SparseMatrix<double> AS( N, N );
      Utility::fill_band( AS, 0, -30.0 );
      Utility::fill_band( AS, -1, 16.0 );
      Utility::fill_band( AS, 1, 16.0 );
      Utility::fill_band( AS, -2, -1.0 );
      Utility::fill_band( AS, 2, -1.0 );

      cout << N << " rows and " << AS.nelts() << " elts. \n";
      DenseVector<double> BS( N, D );

      SparseMatrix<double> aS( AS );
      DenseVector<double> bS( BS );

      #if defined(PETSC_D)
        cout << " Using the PETSc sparse routine:\n";
        SparseLinearSystem<double> sparse_system( &aS, &bS, "petsc" );
      #else
        cout << " Using the native Gauss-Jordan sparse routine:\n";
        SparseLinearSystem<double> sparse_system( &aS, &bS, "native" );
      #endif

      #ifdef TIME
      timer.start();
      do
      {
      #endif
        // reset the matrix
        aS = AS;
        bS = BS;
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
        timer.counter()++;
      }
      while ( timer.get_time() < 5000.0 );
      timer.stop();
      timer.print();
      timer.reset();
      #endif

      bS.sub( bD );

      if ( std::abs( bS.two_norm() ) > tol )
      {
        cout << " \033[1;31;48m * Sparse solver does not give same result as dense solver \033[0m\n";
        failed = true;
      }
      else
      {
        cout << " Sparse solver agrees with the dense solver.\n";
      }

      if ( failed )
      {
        cout << " || dense - sparse ||_2 = " << std::abs( bS.two_norm() ) << "\n";
      }
    } // end sparse banded comparison
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

  PetscFinalize();

}
