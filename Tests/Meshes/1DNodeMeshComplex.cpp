/// \file 1DNodeMeshComplex.cpp
/// \ingroup Tests
/// \ingroup Generic
/// A simple check of the OneD_Node_Mesh container that
/// stores nodal data over a given mesh. This checks that
/// mesh data is interpolated and integrated correctly. We
/// write \f$ \cos(x) \f$ to a uniform mesh, then check the
/// integration routines converge appropriately. We then
/// re-interpolate the data onto a non-uniform mesh and
/// check the integral again.

#include <OneD_Node_Mesh.h>
#include <Utility.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  cout.precision( 10 );
  cout << "\n";
  cout << "=== OneD_Node_Mesh: complex data on a complex path  =\n";
  cout << "\n";


  // Number of points
  std::size_t N(11);
  // parameterisation of a complex path 
  DenseVector<double> paramCoords( Utility::uniform_node_vector( 0.0, 1.0, N ) );
  DenseVector<D_complex> z( N, D_complex(0.0,0.0) );
  // the data stored on the complex path
  OneD_Node_Mesh<D_complex, D_complex> Q( z, 1 );
  // i 
  D_complex eye(0.0,1.0);
  
  // Set the variable values to be defined by a Cosine.
  for ( std::size_t i = 0; i < N; ++i ) {
    double s = paramCoords[i];
    Q.coord(i) = s + eye*s*(1-s);
    // mesh stores 
    Q( i, 0 ) = sin( Q.coord(i) );
  }
  Q.dump_gnu("./data.dat");
  Q.dump();

  OneD_Node_Mesh<D_complex,D_complex> P( "./data.dat", N, 1 );
  P.dump();
  
}
