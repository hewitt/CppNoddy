/// \file EVPRayleigh_lapack.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the Rayleigh problem for values \f$ c \f$
/// that satisfy :
/// \f[ (U_B(y)-c) (\phi''(y) - \alpha^2 \phi(y)) - U_B''(y) \phi(y) = 0\,, \f]
/// subject to \f$ \phi( 0 ) = \phi( 2\pi ) = 0 \f$; it determines the
/// critical wavenumber \f$\alpha \f$ such that \f$ c_i=0 \f$ for \f$ U_B(y)=\sin(y) \f$.
/// The test compares the critical wavenumber with the predicted value of \f$ \sqrt(3)/2 \f$.

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
    // base flow profile
    D_complex U( const D_complex& y )
    {
      return std::sin( y );
    }
    // curvature of the baseflow
    D_complex Udd( const D_complex& y )
    {
      return - std::sin( y );
    }
  }
}

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== EVP: Rayleigh modes, Tollmien's example =========\n";
  cout << "\n";

  Example::alpha = 0.8;    // the wavenumber
  double tol = 1.e-5;      // tolerance for the test
  double left = 0.0;       // from y = 0
  double right = 2 * M_PI; // to y= 2*pi
  unsigned N( 801 );
  // a real distribution of nodes
  DenseVector<double> r_nodes( Utility::power_node_vector( left, right, N, 1.0 ) );
  DenseVector<D_complex> c_nodes( r_nodes );
  // make a distribution of nodes in the complex plane
  for ( unsigned i = 0; i < N; ++i )
  {
    D_complex y( r_nodes[ i ] );
    c_nodes[ i ] -= .2 * D_complex( 0.0, 1.0 ) * y * std::exp( - y );
  }

  // make a base flow on the complex distribution of nodes
  OneD_Node_Mesh<D_complex, D_complex> base( c_nodes, 2 );
  for ( unsigned i = 0; i < c_nodes.size(); ++i )
  {
    D_complex y = c_nodes[ i ];
    base( i, 0 ) = Example::U( y );
    base( i, 1 ) = Example::Udd( y );
  }

  // make the Rayleigh EVP, with 'CHANNEL' boundary conditions
  HST::Rayleigh<D_complex> my_ray( base, Example::alpha, "CHANNEL" );
  // do a global solve
  my_ray.global_evp();

  // for ( unsigned i = 0; i < my_ray.eigenvalues().size(); ++i )
  // {
  //   if ( my_ray.eigenvalues()[ i ].imag() > -0.05 )
  //   {
  //     cout << i << " " << my_ray.eigenvalues()[ i ] << "\n";
  //   }
  // }

  unsigned i_ev = 378;
  std::cout << "locally solved\n";
  // on my machine the eigenvalue is number 378 -- i'm not sure if this
  // could potentially change on other machines depending on LAPACK. We'll see.
  my_ray.iterate_to_neutral( i_ev );

  if ( std::abs( my_ray.alpha() - .5*sqrt( 3. ) ) < tol )
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
  cout << "\033[1;31;48m  * FAILED \033[0m\n";
  cout << "    Final error in critical wavenumber = " << std::abs( my_ray.alpha() - .5*sqrt( 3. ) ) << "\n";
  return 1;

}
