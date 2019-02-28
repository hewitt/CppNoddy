/// \file SparseVector_petscd.cpp
/// \ingroup Tests
/// \ingroup Vector
/// A superficial sanity check of one_norm and
/// vector arithmetic for the sparse vector class.

#include <SparseVector.h>
#include "../Utils_Fill.h"

using namespace std;
using namespace CppNoddy;

int main()
{

  cout << "\n";
  cout << "=== Vector: sparse class example  ===================\n";
  cout << "\n";

  SparseVector<double> vecA( 100 );
  Utils_Fill::fill_random( vecA, 25 );
  double oneA = vecA.one_norm();

  SparseVector<double> vecB( 100 );
  Utils_Fill::fill_random( vecB, 25 );
  double oneB = vecB.one_norm();

  SparseVector<double> vecC( 100 );
  vecC = vecA * 2. + vecB * 2.;

  double tol = 1.e-13;
  if ( vecC.one_norm() - ( 2 * oneA + 2 * oneB ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }
  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;

}



