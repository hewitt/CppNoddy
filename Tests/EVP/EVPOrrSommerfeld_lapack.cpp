/// \file EVPOrrSommerfeld_lapack.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the following linear eigenvalue problem for values \f$ c \f$
/// that satisfy :
/// \f[ \phi''(y) - \alpha^2 \phi(y) - \psi(y) = 0\,, \f]
/// \f[ \psi''(y) - \alpha^2 \psi(y) - i \alpha Re \left \{ ( U(y) - c ) \psi(y) - U''(y) \phi \right \} = 0\,, \f]
/// subject to \f$ \phi(\pm 1) = \phi'(\pm 1) = 0 \f$ where
/// \f$ \alpha = 1.02 \f$, \f$ Re = 5772.2 \f$ and \f$ U(y) = 1 - y^2 \f$.
/// The matrix problem is constructed manually in this case, using second-order
/// finite differences.
/// These values approximately correspond to the first neutral temporal mode
/// in plane Poiseuille flow, therefore the test to be satisfied is that an eigenvalue
/// exists with \f$ \Im ( c ) \approx 0 \f$.

#include <EVP_bundle.h>

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Temporal spectra of the Orr-Sommerfeld eqn =\n";
  cout << "===  with a matrix problem assembled by hand.\n";
  cout << "\n";


  // discretise with these many nodal points
  const std::size_t nodes( 301 );
  // we'll solve as TWO second order problems
  const std::size_t N( 2 * nodes );
  // domain boundaries
  const double left = -1.0;
  const double right = 1.0;
  // spatial step for a uniform mesh
  const double d = ( right - left ) / ( nodes - 1 );

  // matrices for the EVP, initialised with zeroes
  DenseMatrix<D_complex> a( N, N, 0.0 );
  DenseMatrix<D_complex> b( N, N, 0.0 );

  // streamwise wavenumber and Reynolds number
  const double alpha ( 1.02 );
  const double Re ( 5772.2 );
  const D_complex I( 0.0, 1.0 );

  // boundary conditions at the left boundary
  a( 0, 0 ) = 1.0;           // phi( left ) = 0
  a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
  a( 1, 2 ) = 2.0 / d;
  a( 1, 4 ) = -0.5 / d;
  // fill the interior nodes
  for ( std::size_t i = 1; i <= nodes - 2; ++i )
  {
    // position in the channel
    const double y = left + i * d;
    // Poiseuille flow profile
    const double U = ( 1.0 - y * y );
    const double Udd = -2.0;

    // the first quation at the i'th nodal point
    std::size_t row = 2 * i;
    a( row, row ) = -2.0 / ( d * d ) - alpha * alpha;
    a( row, row - 2 ) = 1.0 / ( d * d );
    a( row, row + 2 ) = 1.0 / ( d * d );
    a( row, row + 1 ) = -1.0;

    row += 1;
    // the second equation at the i'th nodal point
    a( row, row ) = -2.0 / ( d * d ) - alpha * alpha - I * alpha * Re * U;
    a( row, row - 2 ) = 1.0 / ( d * d );
    a( row, row + 2 ) = 1.0 / ( d * d );
    a( row, row - 1 ) = I * alpha * Re * Udd;

    b( row, row ) = - I * alpha * Re;
  }
  // boundary conditions at right boundary
  a( N - 2, N - 2 ) = 1.5 / d;
  a( N - 2, N - 4 ) = -2.0 / d;
  a( N - 2, N - 6 ) = 0.5 / d; // psi'( right ) = 0
  a( N - 1, N - 2 ) = 1.0;     // psi( right ) = 0
  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  DenseLinearEigenSystem<D_complex> system( &a, &b );

  try
  {
    system.eigensolve();
  }
  catch (const std::runtime_error &error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }
  // tag any eigenvalues with imaginary part > -0.1
  system.set_shift( D_complex( 0.0, -0.1 ) );
  system.tag_eigenvalues_upper( + 1 );
  lambdas = system.get_tagged_eigenvalues();
  //lambdas.dump();
  double min_growth_rate( lambdas[ 0 ].imag() );
  // make sure we have a near neutral mode
  const double tol = 1.e-3;

  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  TrackerFile spectrum( "./DATA/spectrum.dat" );
  spectrum.push_ptr( &lambdas, "evs" );
  spectrum.update();
  if ( std::abs( min_growth_rate ) < tol )
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

  cout << "\033[1;31;48m  * FAILED \033[0m\n";
  cout << "    Final error = " << min_growth_rate << "\n";
  return 1;
}
