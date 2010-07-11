/// \file 1DNodeMesh.cpp
/// \ingroup Examples
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
  cout << "=== OneD_Node_Mesh: wrapper =========================\n";
  cout << "\n";
  cout << " The test checks the convergence rate for an integral \n";
  cout << " over an increasing number of mesh points. \n\n";

  double error_trapezium = 1.0;
  double error_Simpson = 1.0;
  double error_remesh = 1.0;
  double error_diffroot = 0.0;

  {
    double old_integral = 0.0;
    double integral = 0.0;
    std::size_t n = 2;
    const unsigned jmax = 8;
    for ( unsigned j = 1; j <= jmax ; ++j )
    {
      // A mesh of n points
      n *= 2;
      n += 1;
      // make a uniform mesh
      OneD_Node_Mesh<double> Q( Utility::uniform_node_vector( 0.0, 1.0, n ), 1 );

      // Set the variable values to be defined by a Cosine.
      for ( std::size_t i = 0; i < n; ++i )
      {
        // mesh stores cos(x)
        Q( i, 0 ) = cos( Q.coord( i ) );
      }

      // cos integrates to sin(x) & limits are 0 to 1
      integral = std::abs( Q.integral2( 0 ) - sin( 1.0 ) );
      if ( j > 1 )
      {
        // Check the lin integral method.
        cout << " n = " << n << " |Integral Error| = " << integral
             << " Ratio = " << old_integral / integral << "\n";
      }
      old_integral = integral;
    }
    error_trapezium = integral;
  }

  cout << "\n Checking Simpson integration routine \n\n";

  {
    double old_integral = 0.0;
    double integral = 0.0;
    std::size_t n = 2;
    const unsigned jmax = 8;
    for ( unsigned j = 1; j <= jmax ; ++j )
    {
      // A mesh of n points
      n *= 2;
      n += 1;
      // make a uniform mesh
      OneD_Node_Mesh<double> Q( Utility::uniform_node_vector( 0.0, 1.0, n ), 1 );

      // Set the nodal values to be defined by a Cosine.
      for ( std::size_t i = 0; i < n; ++i )
      {
        // mesh contains cos(x)
        Q( i, 0 ) = cos( Q.coord( i ) );
      }

      // mesh integrates to give sin(x) with limits 0 to 1
      integral = std::abs( Q.integral4( 0 ) - sin( 1.0 ) );
      if ( j > 1 )
      {
        // Check the lin integral method.
        cout << " n = " << n << " |Integral Error| = " << integral
             << " Ratio = " << old_integral / integral << "\n";
      }
      old_integral = integral;
    }
    error_Simpson = integral;
  }

  std::cout << "\n Checking the mesh linear interpolation. \n";
  {
    const std::size_t n = 501;

    // make a nonuniform mesh distribution
    OneD_Node_Mesh<D_complex> Q( Utility::power_node_vector( 0.0, 1.0, n, 2.0 ), 1 );

    for ( std::size_t i = 0; i < n; ++i )
    {
      // set mesh to contain cos(x) but this time with non-uniform x nodes
      Q( i, 0 ) = std::cos( Q.coord( i ) );
    }

    // integrate over the nonuniform mesh
    const double Inonuniform = std::abs( Q.integral2( 0 ) );

    // make a uniform mesh of nodes
    DenseVector<double> Xuniform = Utility::uniform_node_vector( 0.0, 1.0, n );

    // remesh the nonlinear distribution to a uniform mesh
    Q.remesh1( Xuniform );

    const double Iuniform = std::abs( Q.integral2( 0 ) );
    error_remesh = std::abs( Inonuniform - Iuniform );

    std::cout << " |error| = " << error_remesh << "\n";
  }


  std::cout << "\n Checking mesh root interpolation. \n";
  {
    const std::size_t n = 401;

    OneD_Node_Mesh<double> F( Utility::uniform_node_vector( 0.0, 10.0, n ), 1 );

    for ( std::size_t i = 0; i < n; ++i )
    {
      // put sin(x) into var 0
      F( i, 0 ) = std::sin( F.coord( i ) );
    }


    // find the roots of sin(x) with x in 0 to 10
    {
      DenseVector<double> roots( F.find_roots1( 0 ) );
      for ( std::size_t i = 0; i < roots.size(); ++i )
      {
        // roots should be ordered, so just check that they are at n * pi
        error_diffroot = std::max( std::abs( roots[ i ] - ( i + 1 ) * M_PI ) , error_diffroot );
      }
    }

  }

  if ( ( error_trapezium > 1.e-6 ) || ( error_Simpson > 1.e-12 )
       || ( error_remesh > 1.e-6 ) || ( error_diffroot > 1.e-6 ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
