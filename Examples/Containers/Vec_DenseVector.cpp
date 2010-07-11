/// \file Vec_DenseVector.cpp
/// \ingroup Examples
/// \ingroup Containers
/// Some simple sanity checks for the NVector class
/// with both double and complex types.

#include <algorithm>

#include <Types.h>
#include <Timer.h>
#include <Utility.h>
#include <Functors.h>

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Vector: (dense) double/complex Example  ============\n";
  cout << "\n";

  const unsigned N = 1000000; // size of vectors

  cout << " Using vectors of size  " << N << " \n";

  const double tol = 1.e-13;
  bool failed = false;

  Utility::time_seed();
  // a real vector of random entries
  DenseVector<double> V( N, 0.0 );
  Utility::fill_random( V );
  // a complex vector initialised with the same real vector
  DenseVector<D_complex> CExample( V );

  // another real vector of random entries
  DenseVector<double> DExample( N, 0.0 );
  Utility::fill_random( DExample );

#ifdef TIME

  Timer timer;
  const unsigned M = 1000;    // number of repeats
  // time/check BLAS & native dot products
  timer.start();
  for ( unsigned i = 0; i < M; ++i )
  {
    Utility::dot( V, DExample );
  }
  timer.stop();
  std::cout << " dot product done " << M << " times \n";
  timer.print();
#endif

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

#ifdef TIME
  timer.start();
  do
  {
#endif
    CExample.one_norm();
    CExample.two_norm();
    CExample.inf_norm();
#ifdef TIME

    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.stop();
  timer.print();
#endif

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

#ifdef TIME
  timer.start();
  do
  {
#endif

    DExample.one_norm();
    DExample.two_norm();
    DExample.inf_norm();

#ifdef TIME

    timer.counter()++;
  }
  while ( timer.get_time() < 5000.0 );
  timer.stop();
  timer.print();
#endif

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
