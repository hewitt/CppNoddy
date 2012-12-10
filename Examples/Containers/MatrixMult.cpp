/// \file MatrixMult.cpp
/// \ingroup Examples
/// \ingroup Containers
/// Do some matrix multiplications and compare the
/// native N^3 multiply with the BLAS implementation. The
/// test is split between a simple \f$ 2\times 3 \f$,
/// \f$ 3 \times 2\f$ multiply, followed by larger
/// matrices with random contents for a timing test.


#include <Timer.h>
#include <Types.h>
#include <Utility.h>


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Matrix: Native/BLAS multiplication  =============\n";
  cout << "\n";

  DenseMatrix<double> a( 2, 3, 0.0 );
  DenseMatrix<double> b( 3, 2, 0.0 );
  DenseMatrix<double> c( 2, 2, 0.0 );
  a( 0, 0 ) = 1.0;
  a( 0, 1 ) = 2.0;
  a( 0, 2 ) = 3.0;
  a( 1, 0 ) = 4.0;
  a( 1, 1 ) = 5.0;
  a( 1, 2 ) = 6.0;
  b( 0, 0 ) = 1.0;
  b( 0, 1 ) = 2.0;
  b( 1, 0 ) = 3.0;
  b( 1, 1 ) = 4.0;
  b( 2, 0 ) = 5.0;
  b( 2, 1 ) = 6.0;
  DenseMatrix<double> answer( 2, 2, 0.0 );
  answer( 0, 0 ) = 22.0;
  answer( 0, 1 ) = 28.0;
  answer( 1, 0 ) = 49.0;
  answer( 1, 1 ) = 64.0;

  c = a.multiply( b );

  const double tol = 1.e-10;
  bool failed = false;
  c.sub( answer );
  if ( c.inf_norm() > tol )
  {
    std::cout << " Infinity norm of error = " << c.inf_norm() << "\n";
    std::cout << " Native method : Simple (2x3) * (3x2) matrix mult failed\n";
    failed = true;
  }
  else
  {
    std::cout << " Native method : Simple (2x3) * (3x2) matrix mult passed \n";
  }

  // if LAPACK Libs are present, then do the same test
#ifdef LAPACK
  DenseMatrix<double> cblas( 2, 2, 0.0 );
  cblas = Utility::multiply( a, b );
  cblas.sub( answer );
  if ( cblas.inf_norm() > tol )
  {
    std::cout << " Infinity norm of error = " << cblas.inf_norm() << "\n";
    std::cout << " BLAS : Simple (2x3) * (3x2) matrix mult failed \n";
    failed = true;
  }
  else
  {
    std::cout << " BLAS : Simple (2x3) * (3x2) matrix mult test passed \n";
  }
#endif

#ifdef TIME
  Timer timer;
  for ( int N = 128; N <= 2048 ; N *= 2 )
  {
    std::cout << "\n --- Filling arrays of size " << N << "x" << N << "\n";

    // random array data
    DenseMatrix<double> A( N, N, 0.0 );
    Utility::fill_random( A );
    DenseMatrix<double> B( N, N, 0.0 );
    Utility::fill_random( B );
    DenseMatrix<double> C( N, N, 0.0 );

    timer.start();
    C = A.multiply( B );
    timer.stop();
    std::cout << "\n Native N^3 multiplication method :\n";
    timer.print();
    timer.reset();

#ifdef LAPACK

    DenseMatrix<double> Cblas( N, N, 0.0 );
    timer.start();
    Cblas = Utility::multiply( A, B );
    timer.stop();
    std::cout << "\n BLAS multiplication method :\n";
    timer.print();
    timer.reset();
    // Check the difference between the methods
    Cblas.sub( C );
    if ( Cblas.inf_norm() > tol )
    {
      std::cout << " Infinity norm of error = " << Cblas.inf_norm() << "\n";
      std::cout << " BLAS & native matrix multiplication disagree! \n";
      failed = true;
    }
#endif // LAPACK

  }
#endif // TIME

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
