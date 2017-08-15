/// \file Matrix_CompareNative.cpp
/// \ingroup Examples
/// \ingroup Containers
/// A quick check that the overhead associated with the matrix container
/// class is less than 5% compared to a native array.

#include <algorithm>

#include <Types.h>
#include <Timer.h>
#include <Utility.h>
#include <TwoD_Node_Mesh.h>

using namespace CppNoddy;
using namespace std;

int main()
{


  cout << "\n";
  cout << "=== Matrix: compare access speeds to native array ===\n";
  cout << "\n";


  const std::size_t L = 5000;
  const int loops( 100 );
  // DenseMatrix object
  DenseMatrix<double> A( L, L, 0.0 );
  // A TwoD_Node_Mesh object (coordinates of the mesh don't matter)
  const DenseVector<double> coords( Utility::uniform_node_vector( 0.0, 1.0, L ) );
  TwoD_Node_Mesh<double> M( coords, coords, 1 );
  // native array
  double* B;
  B = new double[ L * L ];
  // put some junk into the DenseMatrix
  Utility::fill_random( A );
  // copy the matrix to a native array & the mesh object
  for ( std::size_t row = 0; row < L; ++row )
  {
    for ( std::size_t col = 0; col < L; ++col )
    {
      B[ row * L + col ] = A( row, col );
      M( row, col, 0 ) = A( row, col );
    }
  }

  cout << " DenseMatrix<double> scaling on a per-element basis via access operator.\n";
  Timer timer;
  timer.start();
  for ( int i = 0; i < loops; ++i )
  {
    for ( std::size_t row = 0; row < L; ++row )
    {
      for ( std::size_t col = 0; col < L; ++col )
      {
        A( row, col ) *= 2.0;
      }
    }
  }
  timer.stop();
  double timeA = timer.get_time();
  timer.print();
  timer.reset();

  cout << "\n TwoD_Node_Mesh<double> scaling on a per-element basis via access operator.\n";
  timer.start();
  for ( int i = 0; i < loops; ++i )
  {
    for ( std::size_t row = 0; row < L; ++row )
    {
      for ( std::size_t col = 0; col < L; ++col )
      {
        M( row, col, 0 ) *= 2.0;
      }
    }
  }
  timer.stop();
  double timeM = timer.get_time();
  timer.print();
  timer.reset();

  cout << "\n Native array scaling on a per-element basis.\n";
  timer.start();
  for ( int i = 0; i < loops; ++i )
  {
    for ( std::size_t row = 0; row < L; ++row )
    {
      for ( std::size_t col = 0; col < L; ++col )
      {
        B[ row * L + col ] *= 2.0;
      }
    }
  }
  timer.stop();
  double timeB = timer.get_time();
  timer.print();
  timer.reset();

  delete[] B;  

  bool failed( false );
  // fail if there is more than a 5% overhead between DenseMatrix & native
  if ( ( timeA - timeB ) / ( 0.5 * ( timeA + timeB ) ) > 0.05 )
  {
    failed = true;
    cout << "The % slow-down for a DenseMatrix was " <<
         100.0 * ( timeA - timeB ) / ( 0.5 * ( timeA + timeB ) ) << "\n";
  }
  // fail if there is more than a 5% overhead between DenseMatrix & TwoD_Node_Mesh
  if ( std::abs( timeM - timeB ) / ( 0.5 * ( timeM + timeB ) ) > 0.05 )
  {
    failed = true;
    cout << "The % slow-down for a TwoD_Node_Mesh was " <<
         100.0 * ( timeM - timeB ) / ( 0.5 * ( timeM + timeB ) ) << "\n";
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
