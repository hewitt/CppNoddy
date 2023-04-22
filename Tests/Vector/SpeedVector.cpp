/// \file DenseVector.cpp
/// \ingroup Tests
/// \ingroup Vector
/// Some simple sanity checks for the NVector class
/// with both double and complex types.

#include <algorithm>

#include <Types.h>
#include <Timer.h>
#include <Utility.h>
#include "../Utils_Fill.h"

using namespace CppNoddy;
using namespace std;

int main()
{


  const unsigned N = 250000000; // size of vectors
  cout << " Using vectors of size  " << N << " \n";

  std::vector<double> A( N, 0.0 );
    
  double* B;
  B = new double[ N ];

  DenseVector<double> C( N, 0.0 );
  
  Timer timer;
  cout << "\n Filling a std::vector via [].\n";
  timer.start();
  for ( std::size_t i = 0; i < N; ++i ) {
    A[i] = i;
  }
  timer.stop();
  double timeSTL = timer.get_time();
  timer.print();
  timer.reset();

  cout << "\n Filling a native array via [].\n";
  timer.start();  
  for ( std::size_t i = 0; i < N; ++i ) {
    B[i] = i;
  }
  timer.stop();
  double timeNative = timer.get_time();
  timer.print();
  timer.reset();
  
  
  cout << "\n Filling a DenseVector via [].\n";
  timer.start();  
  for ( std::size_t i = 0; i < N; ++i ) {
    C[i] = i;
  }
  timer.stop();
  double timeDenseVec = timer.get_time();
  timer.print();
  timer.reset();
  
  cout << "\nDenseVector slow down is " << 100*(timeDenseVec-timeNative)/timeNative << "\n";
  cout << "STL:vector slow down is " << 100*(timeSTL-timeNative)/timeNative << "\n";
}
