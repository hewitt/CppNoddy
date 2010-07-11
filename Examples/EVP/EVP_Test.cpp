#include <EVP_bundle.h>
#include <HST.h>

namespace CppNoddy
{
  namespace Example
  {
    // complex base flow in complex plane
		OneD_Node_Mesh<D_complex, D_complex> baseflow;
    // Rayleigh wavenumber
		double alpha;
    double eps;
    double gamma;
    double y0;
    // base flow profile
    double U( const double& y )
    {
      return std::exp( -y ) + eps * std::exp( -gamma*(y-y0)*(y-y0) );
    }   
    // curvature of the baseflow
    double Udd( const double& y )
    {
      return std::exp( -y ) - 2 * gamma * eps * std::exp( -gamma*(y-y0)*(y-y0) ) 
        + 4 * gamma * gamma * (y-y0) * (y-y0) * eps * std::exp( -gamma*(y-y0)*(y-y0) );
    }
  }
}

using namespace CppNoddy;
using namespace std;

int main()
{

  Example::alpha = 0.26;   // the wavenumber
  Example::eps = 0.035;
  Example::gamma = 2;
  Example::y0 = 2;
  double left = 0;      
  double right = 12;       
	unsigned N( 601 );
  // a real distribution of nodes  
	DenseVector<double> r_nodes( Utility::power_node_vector( left, right, N, 1.0 ) );
    
  // make a base flow on the complex distribution of nodes  
  OneD_Node_Mesh<double> base( r_nodes, 2 );
  for ( unsigned i = 0; i < r_nodes.size(); ++i )
  {
    double y = r_nodes[ i ];
    base( i, 0 ) = Example::U( y );
    base( i, 1 ) = Example::Udd( y );
  }

  //do
  //{  
    // make the Rayleigh EVP
    HST::Rayleigh<double> my_ray( base, Example::alpha, "BL" );
    // do a global solve
    my_ray.global_evp();
    
    for ( unsigned i = 0; i < my_ray.eigenvalues().size(); ++i )
    {
      if ( my_ray.eigenvalues()[ i ].imag() > 0.00001 )
      {
        cout << Example::alpha << " " << real(my_ray.eigenvalues()[ i ]) << " " << imag(my_ray.eigenvalues()[ i ]) << "\n";
      } 
    }
    cout << "# " << Example::alpha << "\n";
    //Example::alpha -= 0.01;
  //} while ( Example::alpha > 0.0 );

  HST::Orr_Sommerfeld my_os( base, Example::alpha, 80000.0 );
  my_os.global_evp();
  for ( unsigned i = 0; i < my_os.eigenvalues().size(); ++i )
  {
    if ( my_os.eigenvalues()[ i ].imag() > 0.00001 )
    {
      cout << Example::alpha << " " << real(my_os.eigenvalues()[ i ]) << " " << imag(my_os.eigenvalues()[ i ]) << "\n";
    } 
  }

  TrackerFile temp( "./DATA/evs.dat" );
  temp.push_ptr( &my_os.eigenvalues(), "c" );
  temp.update();
  
}
