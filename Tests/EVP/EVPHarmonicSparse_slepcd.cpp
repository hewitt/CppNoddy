/// \file EVPHarmonicSparse_slepcd.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the harmonic equation
/// \f[ f''(x) + \lambda f(x) = 0 \f]
/// as an eigenvalue problem for \f$ \lambda \f$ over the unit domain
/// with homogeneous boundary conditions for \f$ f(x) \f$, returning
/// any eigenvalue(s) with absolute value less than 10. The resulting
/// eigenvalues should approach \f$ m^2 \pi^2 \f$ as \f$n \to \infty \f$ with
/// \f$ O(\Delta^2)\f$ corrections where \f$ \Delta = 1/(n-1) \f$ and \f$ n \f$ is the number
/// of nodal points in the FD representation; a 2nd order
/// central difference representation of
/// \f[ f''(x_i) + \lambda f(x_i) = \frac{ f_{i-1} - 2f_i + f_{i+1} }{ \Delta^2 } + \lambda f_i + O(\Delta^2) \f]
/// is used. The matrix problem is solved for a subset of eigenvalues within a specified distance
/// of the origin of the comlplex plane by calling the SLEPc generalised eigenvalue routine.

#include <Utility.h>
#include <Timer.h>
#include <SparseLinearEigenSystem.h>
#include <SlepcSession.h>

#include "../Utils_Fill.h"

using namespace CppNoddy;
using namespace std;

int main(int argc, char* argv[])
{
  SlepcSession::getInstance(argc,argv);
  
  cout << "\n";
  cout << "=== EVP: Harmonic equation solved using SLEPc  ======\n";
  cout << "===  with a manually assembled matrix problem.\n";
  cout << "\n";
  
  cout.precision( 12 );
  cout << " Number of nodal points : Leading eigenvalue error : Total CPU time taken (ms) \n";
  bool failed = false;
  size_t N = 4;
  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  Timer timer;
  for ( int i = 2; i < 11; ++i )
  {
    N = ( size_t ) ( std::pow( 2., i ) );
    const double delta = 1. / ( N - 1 );
    const double delta2 = delta * delta;
    // matrix problem
    typedef double PETSC_type;
    SparseMatrix<PETSC_type> a( N, N );
    // Finite difference representation of f''(x)
    // here it's a a tri-diagonal system
    Utils_Fill::fill_band(a, -1, (PETSC_type)(1.0 / delta2));
    Utils_Fill::fill_band(a,  0, (PETSC_type)(-2.0 / delta2));
    Utils_Fill::fill_band(a, 1, (PETSC_type)(1.0 / delta2));
    // overwrite with boundary conditions at f(0) = f(1) = 0
    a( 0, 0 ) = 1.0;
    a( 0, 1 ) = 0.0;
    a( N - 1, N - 1 ) = 1.0;
    a( N - 1, N - 2 ) = 0.0;
    // not a generalised problem - but we'll apply that routine anyway
    // b is the RHS matrix, so it's -I
    SparseMatrix<PETSC_type> b( N, N );
    Utils_Fill::fill_identity( b );
    b.scale( -1.0 );
    b( 0, 0 ) = 0.0;
    b( N - 1, N - 1 ) = 0.0;
    // a vector for the eigenvectors - although we won't use them
    DenseMatrix<D_complex> eigenvectors;

    SparseLinearEigenSystem<PETSC_type> system( &a, &b );
    system.set_calc_eigenvectors( true );
    system.set_nev(4);
    system.set_order( "EPS_TARGET_MAGNITUDE" ); //default

    timer.start();
    try
    {
      system.eigensolve();
    }
    catch (const std::runtime_error &error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }
    system.tag_eigenvalues_disc( + 1, 10. );
    lambdas = system.get_tagged_eigenvalues();
    eigenvectors = system.get_tagged_eigenvectors();
    cout << "    " << N << " : " << lambdas[ 0 ].real() - M_PI * M_PI
         << " : " << timer.get_time() << "\n";
    timer.stop();
  }


  const double tol = 1.e-4;
  if ( abs( lambdas[ 0 ].real() - M_PI * M_PI ) > tol )
    failed = true;

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "    Final error = " << abs( lambdas[ 0 ].real() - M_PI * M_PI ) << "\n";
    return 1;
  }

  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;
}
