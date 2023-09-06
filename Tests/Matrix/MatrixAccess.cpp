/// \file MatrixAccess.cpp
/// \ingroup Tests
/// \ingroup Matrix
/// A quick check that the overhead associated with the matrix container
/// class is less than 5% compared to a native array.

#include <algorithm>

#include <Timer.h>
#include <Types.h>
#include <Utility.h>
#include <TwoD_Node_Mesh.h>

#include "../Utils_Fill.h"

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Matrix: compare access speeds to native array ===\n";
  cout << "\n";

  const std::size_t L = 10000;
  // DenseMatrix object
  DenseMatrix<double> A( L, L, 0.0 );
  // A TwoD_Node_Mesh object (coordinates of the mesh don't matter)
  const DenseVector<double> coords( Utility::uniform_node_vector( 0.0, 1.0, L ) );
  TwoD_Node_Mesh<double> M( coords, coords, 1 );
  // native array
  double* B;
  B = new double[ L * L ];
  // put some junk into the DenseMatrix
  Utils_Fill::fill_random( A );


  // copy the random matrix to a native array & the mesh object
  for ( std::size_t row = 0; row < L; ++row ) {
    for ( std::size_t col = 0; col < L; ++col ) {
      B[ row * L + col ] = A( row, col );
      M( row, col, 0 ) = A( row, col );
    }
  }

  cout << " DenseMatrix<double> ad-hoc op on a per-element basis via access operator.\n";
  Timer timer;
  timer.start();
  for ( std::size_t row = 0; row < L; ++row ) {
    for ( std::size_t col = 0; col < L; ++col ) {
      A( row, col ) = sin( A(row,col) );
    }
  }
  timer.stop();
  double timeA = timer.get_time();
  timer.print();
  timer.reset();

  cout << "\n TwoD_Node_Mesh<double> ad-hoc op on a per-element basis via access operator.\n";
  timer.start();
  for ( std::size_t row = 0; row < L; ++row ) {
    for ( std::size_t col = 0; col < L; ++col ) {
      M( row, col, 0 ) = sin( M(row,col,0) );
    }
  }
  timer.stop();
  double timeM = timer.get_time();
  timer.print();
  timer.reset();

  cout << "\n Native array scaling on a per-element basis.\n";
  timer.start();
  for ( std::size_t row = 0; row < L; ++row ) {
    for ( std::size_t col = 0; col < L; ++col ) {
      B[ row * L + col ] = sin( B[row*L+col] );
    }
  }
  timer.stop();
  double timeB = timer.get_time();
  timer.print();
  timer.reset();

  delete[] B;

  //cout << "The % slow-down for a DenseMatrix was " <<
  //  100.0 * ( timeA - timeB ) / ( 0.5 * ( timeA + timeB ) ) << "\n";

  //cout << "The % slow-down for a TwoD_Node_Mesh was " <<
  //  100.0 * ( timeM - timeB ) / ( 0.5 * ( timeM + timeB ) ) << "\n";
  
  bool failed( false );
  // fail if there is more than a 1% overhead between DenseMatrix & native
  if ( ( timeA - timeB ) / ( 0.5 * ( timeA + timeB ) ) > 0.01 ) {
    failed = true;
    cout << "The % slow-down for a DenseMatrix was " <<
         100.0 * ( timeA - timeB ) / ( 0.5 * ( timeA + timeB ) ) << "%\n";
  }
  // fail if there is more than a 10% overhead between DenseMatrix & TwoD_Node_Mesh
  if ( ( timeM - timeB ) / ( 0.5 * ( timeM + timeB ) ) > 0.10 ) {
    failed = true;      
    cout << "The % slow-down for a TwoD_Node_Mesh was " <<
         100.0 * ( timeM - timeB ) / ( 0.5 * ( timeM + timeB ) ) << "\n";
  }

  if ( failed ) {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  } else {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

}
