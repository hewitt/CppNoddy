/// \file 1DNodeMeshAiry_lapack.cpp
/// \ingroup Tests
/// \ingroup Generic
/// Solves the Airy equation
/// \f[ f''(x) - xf(x) = 0 \f]
/// on the negative real axis over the
/// range [-10,0] then pushes the data into a OneD_GenMesh object
/// and integrates the result. The integral is checked against the
/// table in Abramowitz & Stegun:
/// \f[ \int_{-10}^0 \mbox{Ai}(s) \, \mbox{d}s \approx 0.7656984 \f]
/// The finite-difference solution of
/// the Airy equation and the trapezoidal scheme used in the
/// integration are both second order.

#include <Utility.h>
#include <OneD_Node_Mesh.h>
#include <BandedLinearSystem.h>

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== OneD_Node_Mesh & BandedMatrix: Airy function ====\n";
  cout << "\n";

  size_t n = 8001;                 // number of points
  double l = -10.0;                // domain size
  double d = abs( l ) / ( n - 1 ); // Mesh step
  // a uniform mesh to store the result in
  OneD_Node_Mesh<double> soln( Utility::uniform_node_vector( l, 0.0, n ), 1 );

  BandedMatrix<double> a( n, 3, 0.0 );     // Tridiagonal banded matrix
  DenseVector<double> b( n, 0.0 );         // RHS vector

  a( 0, 0 ) = 1.0;               // BC is that f(-10) = Ai(-10)
  b[ 0 ] = 0.04024123849;        // value from Abramowitz & Stegun
  for ( size_t i = 1; i < n - 1; ++i )
  {
    // set up the f''(x) operator
    a( i, i-1 ) = 1.0/(d*d);
    a( i, i ) = -2.0/(d*d);
    a( i, i+1 ) = 1.0/(d*d);
    // add the -xf(x) term
    a( i, i ) -= soln.coord( i );
  }
  a( n - 1, n - 1 ) = 1.0;       // BC is that f(0) = Ai(0)
  b[ n - 1 ] = 0.3550280539;     // value from Abramowitz & Stegun

  // solve the FD equations
  BandedLinearSystem<double> system( &a, &b, "lapack" );

  try
  {
    system.solve();
  }
  catch (const std::runtime_error &error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }

  // put the result into the solution mesh
  soln.set_vars_from_vector( b );

  // check the integral against Abramowitz & Stegun data
  const double tol = 1.e-5;
  if ( abs( soln.integral2() - 0.7656984 ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout.precision( 10 );
    cout << n << " " << soln.integral2() - 0.7656984 << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

}
