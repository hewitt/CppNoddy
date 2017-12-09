/// \file MatrixSparseSolves_petscd.cpp
/// \ingroup Test
/// \ingroup Matrix
/// Example of the simple linear solvers implemented
/// for sparse matrix objects. A simple
/// \f$ 2 \times 2 \f$ matrix problem is solved using
/// the PETSC_D/Z compiler definitions.


#include <cassert>

#include <Timer.h>
#include <Types.h>
#include <Utility.h>
#include <SparseLinearSystem.h>
#include <petsc.h>
#include <mpi.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  PetscInitialize(NULL,NULL,(char*)0,(char*)0);

  cout << "\n";
  cout << "=== Matrix: Example linear (double) sparse solver  ===\n";
  cout << "\n";

  bool failed = false;
  // tolerance for the test
  const double tol = 1.e-10;

  //
  // SOLVE A SMALL "Sparse"(!) 2X2 REAL SYSTEM
  cout << "=== Matrix: double  ===============================\n";
  //
  SparseMatrix<double> A( 2, 2 );
  DenseVector<double> B( 2, 0.0 );
  A( 0, 0 ) = 1.;
  A( 0, 1 ) = 2.;
  A( 1, 0 ) = 3.;
  A( 1, 1 ) = 4.;
  B[ 0 ] = 5.;
  B[ 1 ] = 11.;

  SparseLinearSystem<double> small_system( &A, &B, "petsc" );

  try
  {
    small_system.factorise();
    small_system.solve_using_factorisation();
  }
  catch ( std::runtime_error )
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
    std::cout << "\033[1;31;48m Simple 2x2 double sparse system was not solved correctly\033[0m\n";
    std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
    failed = true;
  }
  else
  {
    std::cout << " Simple 2x2 double sparse solve works.\n";
  }

  // reset B to be double the previous case
  B[ 0 ] = 10.;
  B[ 1 ] = 22.;
  try
  {
    small_system.solve_using_factorisation();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }
  // double the RHS and double the solution
  answer[ 0 ] = 2.0;
  answer[ 1 ] = 4.0;
  B.sub( answer );
  if ( B.inf_norm() > tol )
  {
    std::cout << "\033[1;31;48m Simple 2x2 double sparse system was not solved correctly\033[0m\n";
    std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
    failed = true;
  }
  else
  {
    std::cout << " Simple 2x2 double sparse solve_using_factorisation works.\n";
  }

  PetscFinalize();

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
