/// \file DistributedVector_petscd.cpp
/// \ingroup Tests
/// \ingroup Vector
/// A superficial sanity check of copy and two_norm

#include <DistributedVector.h>
#include <Types.h>
#include "../Utils_Fill.h"
#include <mpi.h>

using namespace std;
using namespace CppNoddy;


int main(int argc, char *argv[])
{
  PetscSession::getInstance(argc,argv);

  PetscPrintf(PETSC_COMM_WORLD, "\n=== Vector: A distributed (double) example ==========\n\n");

  DistributedVector<double> vecA( 10 );
  
  DenseVector<int> indices( 5, 0 );
  DenseVector<double> values( 5, 1.0);
  for ( auto i = 0; i < 5; i++ ) {
    indices[i] = 2*i;
    values[i] = 2*i;
  }
  vecA.set( indices, values );
  vecA.final_assembly();
  
  auto vecB = vecA;

  vecA.view();
  //vecB.view();
  PetscPrintf(PETSC_COMM_WORLD, "%f", vecA.two_norm() - vecB.two_norm() );
  
  if ( vecA.two_norm() - vecB.two_norm() > 1.e-12 )
  {
    PetscPrintf(PETSC_COMM_WORLD, "\033[1;31;48m  * FAILED \033[0m\n");
    return 1;
  }
  PetscPrintf(PETSC_COMM_WORLD, "\033[1;32;48m  * PASSED \033[0m\n");
  return 0;

}



