/// \file MatrixSolves.cpp
/// \ingroup Examples
/// \ingroup Containers
/// Example of the simple linear solvers implemented
/// for dense, banded and sparse matrix objects. To begin
/// with a simple \f$ 2 \times 2 \f$ matrix problem is solved.
/// Then a penta-diagonal problem is solved using the
/// dense, banded and sparse containers. The native linear Gaussian
/// elimination solvers are used unless the LAPACK/SUPERLU/MUMPS_SEQ compiler
/// options are used, in which case the linear solver phase
/// calls the LAPACK/SUPERLU/MUMPS_SEQ library as appropriate.

#include <cassert>

#include <Timer.h>
#include <Types.h>
#include <Utility.h>
#include <DenseLinearSystem.h>
#include <BandedLinearSystem.h>
#include <SparseLinearSystem.h>

#ifdef MUMPS_SEQ
#include <MPIinit.h>
#endif

using namespace CppNoddy;
using namespace std;

int main()
{

#ifndef MUMPS_SEQ
  cout << "\n";
  cout << "=== Matrix: Example MUMPS MPI linear solver  ========\n";
  cout << "\n";
  cout << " MUMPS/MPI support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
#else
  int myid, size;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  if ( myid == 0 )
  {
    cout << "\n";
    cout << "=== Matrix: Example MUMPS MPI linear solver  ======\n";
    cout << "\n";
  }
  MPI_Barrier( MPI_COMM_WORLD );
  cout << " MPI instance " << myid << " of " << size << " is running.\n";

  bool failed = false;
  // tolerance for the test
  const double tol = 1.e-10;

  const unsigned N = 511;
  const double D = 12 * ( 1. / ( N - 1 ) ) * ( 1. / ( N - 1 ) );
  // SOLVE the BANDED REAL SYSTEM as above but as a sparse system using native and superlu

  DenseMatrix<double> AD( N, N, 0.0 );
  DenseVector<double> BD( N, D );
  Utility::fill_band( AD, 0, -30.0 );
  Utility::fill_band( AD, -1, 16.0 );
  Utility::fill_band( AD, 1, 16.0 );
  Utility::fill_band( AD, -2, -1.0 );
  Utility::fill_band( AD, 2, -1.0 );
  DenseLinearSystem<double> dense_system( &AD, &BD, "lapack" );
  dense_system.solve();

  //
  // SOLVE the BANDED REAL SYSTEM as above but as a sparse system using mumps_seq
  //
  if ( myid == 0 )
  {
    cout << " Comparing the sparse matrix solver solution to dense solver: ";
  }
  SparseMatrix<double> AS( N, N );
  Utility::fill_band( AS, 0, -30.0 );
  Utility::fill_band( AS, -1, 16.0 );
  Utility::fill_band( AS, 1, 16.0 );
  Utility::fill_band( AS, -2, -1.0 );
  Utility::fill_band( AS, 2, -1.0 );
  if ( myid == 0 )
  {
    cout << N << " rows and " << AS.nelts() << " elts. \n";
  }
  DenseVector<double> BS( N, D );

  if ( myid == 0 )
  {
    cout << " Using the MUMPS_SEQ sparse routine:\n";
  }
  SparseLinearSystem<double> sparse_system( &AS, &BS, "mumps_seq" );

  sparse_system.solve();
  BS.sub( BD );

  if ( myid == 0 )
  {
    if ( std::abs( BS.two_norm() ) > tol )
    {
      cout << " \033[1;31;48m * Sparse solver does not give same result as dense solver \033[0m\n";
      failed = true;
    }
    else
    {
      cout << " Sparse (mumps_seq) agrees with the dense solver.\n";
    }
    if ( failed )
    {
      cout << " || dense - sparse ||_2 = " << std::abs( BS.two_norm() ) << "\n";
      cout << "\033[1;31;48m  * FAILED \033[0m\n";
    }
    else
    {
      cout << "\033[1;32;48m  * PASSED \033[0m\n";
    }
  }

#endif //check for MUMPS
}
