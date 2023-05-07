/// \file EVPHarmonicComplex_lapack.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the harmonic equation
/// \f[ f''(z) + \lambda f(z) = 0 \f]
/// as an eigenvalue problem for \f$ \lambda \f$ over a path in the
/// complex plane with homogeneous boundary conditions for \f$ f(z) \f$, returning
/// any eigenvalue(s) with absolute value less than 10. The resulting
/// eigenvalues should approach \f$ m^2 \pi^2 \f$ as \f$n \to \infty \f$ with
/// \f$ O(\Delta^2)\f$ corrections where \f$ \Delta = 1/(n-1) \f$ and \f$ n \f$ is the number
/// of nodal points in the FD representation.
/// The complex path is parametrised by a real parameter "x" and a 2nd order
/// central difference representation is used.
/// The matrix problem is solved for all eigenvalues within a specified distance
/// of the origin of the comlplex plane by calling the LAPACK generalised eigenvalue routine.

#include <Utility.h>
#include <EVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    const D_complex eye(0.,1.0);
    // how far to deform the path into the complex plane
    const double delta(0.5);

    // functions that define the complex path z=z(x)
    D_complex z( double x ) {
      return x + eye*delta*x*(1-x);
    }
    // derivative of the path w.r.t. parameter x
    D_complex zx( double x ) {
      return 1.0 + eye*delta*(1-2.*x);
    }
    // 2nd derivative of the path w.r.t. parameter x
    D_complex zxx( double x ) {
      return -eye*2.*delta;
    }
  }
}

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Harmonic equation solved using LAPACK  =====\n";
  cout << "===  with a manually assembled matrix problem.\n";
  cout << "===  The problem is solved along a path in the complex\n";
  cout << "===  plane, deformed away from the real axis.\n";
  cout << "\n";

  cout.precision( 12 );
  cout << " Number of nodal points : Leading eigenvalue error : Total CPU time taken (ms) \n";
  bool failed = false;
  size_t N = 4;
  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  DenseMatrix<D_complex> eigenvectors;
  
  Timer timer;
  for ( int i = 2; i < 11; ++i )
  {
    N = ( size_t ) ( std::pow( 2., i ) );
    const double delta = 1. / ( N - 1 );
    const double delta2 = delta * delta;
    // matrix problem
    DenseMatrix<D_complex> a( N, N, 0.0 );
    DenseMatrix<D_complex> b( N, N, 0.0 );
    // boundary conditions at f(0) = 0
    a( 0, 0 ) = 1.0;
    a( 0, 1 ) = 0.0;
    for ( unsigned j = 1; j < N-1; ++j ) {
      // Finite difference representation of f''(x)
      double x = j*delta;
      a( j, j-1 ) =  (1.0/pow(Example::zx(x),2.0))/delta2 + (Example::zxx(x)/pow(Example::zx(x),3.0))/(2*delta);
      a( j, j)    = -(2.0/pow(Example::zx(x),2.0))/delta2;
      a( j, j+1 ) =  (1.0/pow(Example::zx(x),2.0))/delta2 - (Example::zxx(x)/pow(Example::zx(x),3.0))/(2*delta);

      // not a generalised problem - but we'll apply that routine anyway
      b( j, j )   =  -1.0;
    }
    // boundary conditions at f(1) = 0
    a( N - 1, N - 1 ) = 1.0;
    a( N - 1, N - 2 ) = 0.0;
    // a vector for the eigenvectors - although we won't use them
    DenseLinearEigenSystem<D_complex> system( &a, &b );
    system.set_calc_eigenvectors( true );

    timer.start();
    try
    {
      system.eigensolve();
    } catch (const std::runtime_error &error ) {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }
    
    system.tag_eigenvalues_disc( + 1, 10. );
    lambdas = system.get_tagged_eigenvalues();
    
    eigenvectors = system.get_tagged_eigenvectors();

    cout << "    " << N << " : " << lambdas[ 0 ].real() - M_PI * M_PI << " +i " << lambdas[0].imag() 
         << " : " << timer.get_time() << "\n";
/    timer.stop();
  }

  // real parameterisation of complex path
  DenseVector<double> xNodes( Utility::uniform_node_vector( 0.0, 1.0, N ) );
  // complex path
  DenseVector<D_complex> zNodes( N, D_complex(0.0,0.0) );
  // mesh of complex data along a complex path
  OneD_Node_Mesh<D_complex, D_complex> mesh( zNodes, 1 );
  int index(0); // first and only eigenvalue returned
  for ( unsigned j = 0; j < N; ++j ) {
    // store the nodes in the complex plane
    mesh.coord(j) = Example::z(xNodes[j]);
    mesh(j,0) = eigenvectors(index,j+0); //f_j
  }
  // write the complex solution on the complex path
  mesh.dump_gnu("/DATA/complexMesh.dat"); 
  const double tol = 1.e-4;
  if ( abs( lambdas[ 0 ].real() - M_PI * M_PI ) > tol )
    failed = true;

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "    Final error = " << abs( lambdas[ 0 ].real() - M_PI * M_PI ) << "\n";
    return 1;
  }

  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;
}


#include <Utility.h>
#include <EVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    const D_complex eye(0.,1.0);
    // how far to deform the path into the complex plane
    const double delta(0.5);

    // functions that define the complex path z=z(x)
    D_complex z( double x ) {
      return x + eye*delta*x*(1-x);
    }
    // derivative of the path w.r.t. parameter x
    D_complex zx( double x ) {
      return 1.0 + eye*delta*(1-2.*x);
    }
    // 2nd derivative of the path w.r.t. parameter x
    D_complex zxx( double x ) {
      return -eye*2.*delta;
    }
  }
}

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Harmonic equation solved using LAPACK  =====\n";
  cout << "===  with a manually assembled matrix problem.\n";
  cout << "===  The problem is solved along a path in the complex\n";
  cout << "===  plane, deformed away from the real axis.\n";
  cout << "\n";

  cout.precision( 12 );
  cout << " Number of nodal points : Leading eigenvalue error : Total CPU time taken (ms) \n";
  bool failed = false;
  size_t N = 4;
  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  DenseMatrix<D_complex> eigenvectors;
  
  Timer timer;
  for ( int i = 2; i < 11; ++i )
  {
    N = ( size_t ) ( std::pow( 2., i ) );
    const double delta = 1. / ( N - 1 );
    const double delta2 = delta * delta;
    // matrix problem
    DenseMatrix<D_complex> a( N, N, 0.0 );
    DenseMatrix<D_complex> b( N, N, 0.0 );
    // boundary conditions at f(0) = 0
    a( 0, 0 ) = 1.0;
    a( 0, 1 ) = 0.0;
    for ( unsigned j = 1; j < N-1; ++j ) {
      // Finite difference representation of f''(x)
      double x = j*delta;
      a( j, j-1 ) =  (1.0/pow(Example::zx(x),2.0))/delta2 + (Example::zxx(x)/pow(Example::zx(x),3.0))/(2*delta);
      a( j, j)    = -(2.0/pow(Example::zx(x),2.0))/delta2;
      a( j, j+1 ) =  (1.0/pow(Example::zx(x),2.0))/delta2 - (Example::zxx(x)/pow(Example::zx(x),3.0))/(2*delta);

      // not a generalised problem - but we'll apply that routine anyway
      b( j, j )   =  -1.0;
    }
    // boundary conditions at f(1) = 0
    a( N - 1, N - 1 ) = 1.0;
    a( N - 1, N - 2 ) = 0.0;
    // a vector for the eigenvectors - although we won't use them
    DenseLinearEigenSystem<D_complex> system( &a, &b );
    system.set_calc_eigenvectors( true );

    timer.start();
    try
    {
      system.eigensolve();
    } catch (const std::runtime_error &error ) {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }
    
    system.tag_eigenvalues_disc( + 1, 10. );
    lambdas = system.get_tagged_eigenvalues();
    
    eigenvectors = system.get_tagged_eigenvectors();

    cout << "    " << N << " : " << lambdas[ 0 ].real() - M_PI * M_PI << " +i " << lambdas[0].imag() 
         << " : " << timer.get_time() << "\n";
/    timer.stop();
  }

  // real parameterisation of complex path
  DenseVector<double> xNodes( Utility::uniform_node_vector( 0.0, 1.0, N ) );
  // complex path
  DenseVector<D_complex> zNodes( N, D_complex(0.0,0.0) );
  // mesh of complex data along a complex path
  OneD_Node_Mesh<D_complex, D_complex> mesh( zNodes, 1 );
  int index(0); // first and only eigenvalue returned
  for ( unsigned j = 0; j < N; ++j ) {
    // store the nodes in the complex plane
    mesh.coord(j) = Example::z(xNodes[j]);
    mesh(j,0) = eigenvectors(index,j+0); //f_j
  }
  // write the complex solution on the complex path
  mesh.dump_gnu("/DATA/complexMesh.dat"); 
  const double tol = 1.e-4;
  if ( abs( lambdas[ 0 ].real() - M_PI * M_PI ) > tol )
    failed = true;

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "    Final error = " << abs( lambdas[ 0 ].real() - M_PI * M_PI ) << "\n";
    return 1;
  }

  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;
}


#include <Utility.h>
#include <EVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    const D_complex eye(0.,1.0);
    // how far to deform the path into the complex plane
    const double delta(0.5);

    // functions that define the complex path z=z(x)
    D_complex z( double x ) {
      return x + eye*delta*x*(1-x);
    }
    // derivative of the path w.r.t. parameter x
    D_complex zx( double x ) {
      return 1.0 + eye*delta*(1-2.*x);
    }
    // 2nd derivative of the path w.r.t. parameter x
    D_complex zxx( double x ) {
      return -eye*2.*delta;
    }
  }
}

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: Harmonic equation solved using LAPACK  =====\n";
  cout << "===  with a manually assembled matrix problem.\n";
  cout << "===  The problem is solved along a path in the complex\n";
  cout << "===  plane, deformed away from the real axis.\n";
  cout << "\n";

  cout.precision( 12 );
  cout << " Number of nodal points : Leading eigenvalue error : Total CPU time taken (ms) \n";
  bool failed = false;
  size_t N = 4;
  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  DenseMatrix<D_complex> eigenvectors;
  
  Timer timer;
  for ( int i = 2; i < 11; ++i )
  {
    N = ( size_t ) ( std::pow( 2., i ) );
    const double delta = 1. / ( N - 1 );
    const double delta2 = delta * delta;
    // matrix problem
    DenseMatrix<D_complex> a( N, N, 0.0 );
    DenseMatrix<D_complex> b( N, N, 0.0 );
    // boundary conditions at f(0) = 0
    a( 0, 0 ) = 1.0;
    a( 0, 1 ) = 0.0;
    for ( unsigned j = 1; j < N-1; ++j ) {
      // Finite difference representation of f''(x)
      double x = j*delta;
      a( j, j-1 ) =  (1.0/pow(Example::zx(x),2.0))/delta2 + (Example::zxx(x)/pow(Example::zx(x),3.0))/(2*delta);
      a( j, j)    = -(2.0/pow(Example::zx(x),2.0))/delta2;
      a( j, j+1 ) =  (1.0/pow(Example::zx(x),2.0))/delta2 - (Example::zxx(x)/pow(Example::zx(x),3.0))/(2*delta);

      // not a generalised problem - but we'll apply that routine anyway
      b( j, j )   =  -1.0;
    }
    // boundary conditions at f(1) = 0
    a( N - 1, N - 1 ) = 1.0;
    a( N - 1, N - 2 ) = 0.0;
    // a vector for the eigenvectors - although we won't use them
    DenseLinearEigenSystem<D_complex> system( &a, &b );
    system.set_calc_eigenvectors( true );

    timer.start();
    try
    {
      system.eigensolve();
    } catch (const std::runtime_error &error ) {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }
    
    system.tag_eigenvalues_disc( + 1, 10. );
    lambdas = system.get_tagged_eigenvalues();
    
    eigenvectors = system.get_tagged_eigenvectors();

    cout << "    " << N << " : " << lambdas[ 0 ].real() - M_PI * M_PI << " +i " << lambdas[0].imag() 
         << " : " << timer.get_time() << "\n";
/    timer.stop();
  }

  // real parameterisation of complex path
  DenseVector<double> xNodes( Utility::uniform_node_vector( 0.0, 1.0, N ) );
  // complex path
  DenseVector<D_complex> zNodes( N, D_complex(0.0,0.0) );
  // mesh of complex data along a complex path
  OneD_Node_Mesh<D_complex, D_complex> mesh( zNodes, 1 );
  int index(0); // first and only eigenvalue returned
  for ( unsigned j = 0; j < N; ++j ) {
    // store the nodes in the complex plane
    mesh.coord(j) = Example::z(xNodes[j]);
    mesh(j,0) = eigenvectors(index,j+0); //f_j
  }
  // write the complex solution on the complex path
  mesh.dump_gnu("/DATA/complexMesh.dat"); 
  const double tol = 1.e-4;
  if ( abs( lambdas[ 0 ].real() - M_PI * M_PI ) > tol )
    failed = true;

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "    Final error = " << abs( lambdas[ 0 ].real() - M_PI * M_PI ) << "\n";
    return 1;
  }

  cout << "\033[1;32;48m  * PASSED \033[0m\n";
  return 0;
}
