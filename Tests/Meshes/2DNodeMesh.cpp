/// \file 2DNodeMesh.cpp
/// \ingroup Tests
/// \ingroup Generic
/// A simple check of the TwoD_Node_Mesh container that
/// stores nodal data over a given mesh. This simply
/// does a remeshing check.

#include <Generic_bundle.h>
#include <TwoD_Node_Mesh.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  cout.precision( 10 );
  cout << "\n";
  cout << "=== TwoD_Node_Mesh: very basic read/write test ======\n";
  cout << "\n";

  std::size_t nx( 11 );
  std::size_t ny( 21 );
  DenseVector<double> x = Utility::uniform_node_vector( -1.0, 1.0, nx );
  DenseVector<double> y = Utility::uniform_node_vector(  0.0, 1.0, ny );
  std::size_t nx2( 31 );
  std::size_t ny2( 41 );
  DenseVector<double> x2 = Utility::uniform_node_vector( -1.0, 1.0, nx2 );
  DenseVector<double> y2 = Utility::uniform_node_vector(  0.0, 1.0, ny2 );


  TwoD_Node_Mesh<double> mesh( x, y, 1 );
  // write
  for ( std::size_t i = 0; i < nx; ++i )
  {
    for ( std::size_t j = 0; j < ny; ++j )
    {
      mesh( i, j, 0 ) = cos( M_PI * x[ i ] ) * sin( M_PI * y[ j ] );
    }
  }

  //mesh.dump_gnu( "./DATA/mesh1.dat" );
  mesh.remesh1( x2, y2 );
  //mesh.dump_gnu( "./DATA/mesh2.dat" );
  mesh.remesh1( x, y );
  //mesh.dump_gnu( "./DATA/mesh3.dat" );

  // read
  for ( std::size_t j = 0; j < ny; ++j )
  {
    for ( std::size_t i = 0; i < nx; ++i )
    {
      mesh( i, j, 0 ) -= cos( M_PI * x[ i ] ) * sin( M_PI * y[ j ] );
    }
  }

  // check diff
  const double diff = mesh.get_var_as_matrix( 0 ).inf_norm();
  if ( diff > 1.e-14 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout.precision( 10 );
    cout << "Difference following remeshing is " << diff << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

}
