/// \file DenseVector.cpp
/// \ingroup Tests
/// \ingroup Vector
/// Some simple sanity checks for the NVector class
/// with both double and complex types.

#include <algorithm>

#include <Types.h>
#include <Timer.h>
#include <Utility.h>
#include <Functors.h>
#include "../Utils_Fill.h"

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Vector: (dense) double/complex Example  =========\n";
  cout << "\n";

  const unsigned N = 1000000; // size of vectors

  cout << " Using vectors of size  " << N << " \n";

  const double tol = 1.e-13;
  bool failed = false;

  Utils_Fill::time_seed();
  // a real vector of random entries
  DenseVector<double> V( N, 0.0 );
  Utils_Fill::fill_random( V );
  // a complex vector initialised with the same real vector
  DenseVector<D_complex> CExample( V );

  // another real vector of random entries
  DenseVector<double> DExample( N, 0.0 );
  Utils_Fill::fill_random( DExample );

  const unsigned M = 100;    // number of repeats
  for ( unsigned i = 0; i < M; ++i )
  {
    Utility::dot( V, DExample );
  }

  cout << " \nComplex vectors \n";
  cout << " Testing norms, nearest_index, max/minabs_index\n";
  if ( abs( CExample.inf_norm() - *min_element( CExample.begin(), CExample.end(), absDiff_predicate<D_complex>( 1.0 ) ) ) >
       tol )
  {
    std::cout << "  - Failed inf_norm/nearest_index Example" << "\n";
    failed = true;
  }
  if ( abs( CExample.inf_norm() - *max_element( CExample.begin(), CExample.end(), abs_predicate<D_complex>() ) )
       > tol )
  {
    std::cout << "  - Failed inf_norm/maxabs_index Example" << "\n";
    failed = true;
  }

  CExample.one_norm();
  CExample.two_norm();
  CExample.inf_norm();

  cout << "\n";
  cout << " Double vectors \n";

  cout << " Testing norms, nearest_index, max/minabs_index\n";
  if ( abs( DExample.inf_norm() - *min_element( DExample.begin(), DExample.end(), absDiff_predicate<double>( 1.0 ) ) )
       > tol )
  {
    std::cout << "  - Failed inf_norm/nearest_index Example" << "\n";
    failed = true;
  }
  if ( abs( DExample.inf_norm() - *max_element( DExample.begin(), DExample.end(), abs_predicate<double>() ) )
       > tol )
  {
    std::cout << "  - Failed inf_norm/maxabs_index Example" << "\n";
    failed = true;
  }

  DExample.one_norm();
  DExample.two_norm();
  DExample.inf_norm();

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
