/// \file MatrixSolves.cpp
/// \ingroup Examples
/// \ingroup Containers
/// Example of the simple linear solvers implemented
/// for sparse matrix objects. A simple
/// \f$ 2 \times 2 \f$ matrix problem is solved using
/// the PETSC_D/Z compiler definitions.


#include <cassert>

#include <Timer.h>
#include <Types.h>
#include <Utility.h>
#include <SparseLinearSystem.h>


using namespace CppNoddy;
using namespace std;

int main()
{
  #if defined(PETSC_Z) || defined(PETSC_D)
  cout << "=== Init of PETSc. ==================================\n";
  PetscInitialize(NULL,NULL,(char*)0,(char*)0);
  #endif

  cout << "\n";
  cout << "=== Matrix: Example linear sparse solver  ===========\n";
  cout << "\n";

  bool failed = false;
  // tolerance for the test
  const double tol = 1.e-10;

#ifdef PETSC_D
  {
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
      assert( false );
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

  }
#endif

#ifdef PETSC_Z
  {
    //
    // SOLVE A SMALL "Sparse"(!) 2X2 COMPLEX SYSTEM
    cout << "=== Matrix: complex  ==============================\n";
    //
    SparseMatrix<D_complex> A( 2, 2 );
    DenseVector<D_complex> B( 2, 0.0 );
    A( 0, 0 ) = 1.;
    A( 0, 1 ) = 2.;
    A( 1, 0 ) = 3.;
    A( 1, 1 ) = 4.;
    B[ 0 ] = 5.;
    B[ 1 ] = 11.;

    SparseLinearSystem<D_complex> small_system( &A, &B, "petsc" );

    try
    {
      // small_system.solve();
      small_system.factorise();
      small_system.solve_using_factorisation();
    }
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      assert( false );
    }
    DenseVector<D_complex> answer( 2, 0.0 );
    answer[ 0 ] = 1.0;
    answer[ 1 ] = 2.0;
    B.sub( answer );
    if ( B.inf_norm() > tol )
    {
      std::cout << "\033[1;31;48m Simple 2x2 complex sparse system was not solved correctly\033[0m\n";
      std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
      failed = true;
    }
    else
    {
      std::cout << " Simple 2x2 complex sparse solve works.\n";
    }

    // double the RHS and try a resolve
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
    answer[ 0 ] = 2.0;
    answer[ 1 ] = 4.0;
    B.sub( answer );
    if ( B.inf_norm() > tol )
    {
      std::cout << "\033[1;31;48m Simple 2x2 complex sparse system was not solved correctly\033[0m\n";
      std::cout << " residual vector's inf_norm = " << B.inf_norm() << "\n";
      failed = true;
    }
    else
    {
      std::cout << " Simple 2x2 complex sparse solve_using_factorisation works.\n";
    }

  }
  cout << "COMPLEX 2x2 TEST ENDED\n";
#endif

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

  #if defined(PETSC_D)||defined(PETSC_Z)
    PetscFinalize();
  #endif
}
