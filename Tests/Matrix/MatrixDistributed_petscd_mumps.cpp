/// \file DistributedMatrixTest_petscz.cpp
/// \ingroup Tests
/// \ingroup Matrix
/// A superficial sanity check

#include <DistributedMatrix.h>
#include <DistributedVector.h>
#include <DistributedLinearSystem.h>
#include <PetscSession.h>
#include <Types.h>
#include <mpi.h>

using namespace std;
using namespace CppNoddy;


int main(int argc, char *argv[])
{
  PetscSession::getInstance(argc,argv);

  PetscPrintf(PETSC_COMM_WORLD, "\n=== Vector: A distributed (double) example ==========\n\n");

  std::size_t Nx(101);
  std::size_t Ny(101);
  std::size_t N(Nx*Ny);
  double deltax = 1.0/(Nx-1);
  double deltaxSq = deltax*deltax;
  double deltay = 1.0/(Ny-1);
  double deltaySq = deltay*deltay;
  DistributedMatrix<double> matA( N, N, 5, 4 );
  DistributedVector<double> vecB( N );

  for(std::size_t j = 0; j < Ny; ++j ) {
    matA.set_elt(j,j,1.0);
    vecB.set_elt(j,0.0);   
  }
  //
  for(std::size_t i = 1; i < Nx-1; ++i ) {
    matA.set_elt( i*Ny, i*Ny, 1.0 );
    vecB.set_elt( i*Ny, 0.0 );   
    for(std::size_t j = 1; j < Ny-1; ++j ) { 
      DenseVector<int> index;
      DenseVector<double> value;
      index.push_back( i*Ny+j-1 );
      value.push_back( 1./deltaySq );
      index.push_back( i*Ny+j );
      value.push_back( -2.0/deltaySq - 2.0/deltaxSq);
      index.push_back( i*Ny+j+1 );
      value.push_back( 1./deltaySq );
      index.push_back( (i-1)*Ny + j );
      value.push_back( 1./deltaxSq );
      index.push_back( (i+1)*Ny + j );
      value.push_back( 1./deltaxSq );      
      matA.set_row( i*Ny+j, index, value );
      vecB.set_elt( i*Ny+j, 1.0 );
    }
    matA.set_elt( i*Ny+Ny-1, i*Ny+Ny-1, 1.0 );
    vecB.set_elt( i*Ny+Ny-1, 0.0 );   
  }
  for(std::size_t j = 0; j < Ny; ++j ) {
    matA.set_elt( Ny*(Nx-1) + j , Ny*(Nx-1) + j, 1.0 );
    vecB.set_elt( Ny*(Nx-1) + j, 0.0 );   
  }
  
  matA.final_assembly();
  vecB.final_assembly();
  
  //matA.view();
  //vecB.view();

  DistributedLinearSystem<double> system( &matA, &vecB );
  system.solve();

  cout << vecB.two_norm() << "\n";
  
  //vecB.view();

}



