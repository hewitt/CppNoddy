/// \file DenseVectorOpOverloading.cpp
/// \ingroup Tests
/// \ingroup Vector
/// Just a quick and simple check that DenseVector operator
/// overloading is functioning.

#include <Types.h>

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Vector: operator overloading check ==============\n";
  cout << "\n";

  DenseVector<double> A1( 3, 1.0 );
  DenseVector<double> B1( 3, 2.0 );
  DenseVector<double> C1( 3, 0.0 );

  const double n1 = 10. / 3.;

  C1 = A1 * 10 + B1 * 5; // C1 has all elts = 20
  C1 *= n1;              // C1 has all elts = 66.66666

  DenseVector<D_complex> A2( 3, 1.0 );
  DenseVector<D_complex> B2( 3, 2.0 );
  DenseVector<D_complex> C2( 3, 0.0 );

  D_complex n2( 10.0, 10.0 );
  D_complex m2( 5.0, 5.0 );
  D_complex l2( 1.0, 1.0 );

  C2 = A2 * n2 + B2 * m2; // C2 has all elts = (20,20)
  C2 *= n1;               // C2 has all elts = (66.6666,66.6666)

  // make a complex version of C1 called C3
  DenseVector<D_complex> C3( C1 );
  C3 = C3 * l2;           // C3 has all elts = (66.6666,66.6666)
  C2 = C2 - C3;

  const double tol = 1.e-14;
  if ( abs( C2.two_norm() ) > tol )
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
