/// \file PoissonCartesian_petscd.cpp
/// \ingroup Tests
/// \ingroup Poisson
/// Solving a Cartesian Poisson problem:
/// \f[ \nabla^2 \psi(x,y) = 2( x^2 + y^2 ) \f]
/// with \f[ \psi(x,\pm 1) = x^2\,, \quad\mbox{and}\quad \psi(\pm 1, y) = y^2 \f]
/// where \f$ (x,y) \in [-1,1]\times[-1,1] \f$.
/// The global problem is solved (one step) and result is compared to the
/// exact solution \f[ \psi(x,y) = x^2y^2 \f].

#include <TwoD_Node_Mesh.h>
#include <SparseLinearSystem.h>
#include <DenseLinearSystem.h>
#include <Utility.h>
#include <PetscSession.h>

namespace CppNoddy
{
  namespace Example
  {
    double source_fn( const std::pair<double,double>& coord ) {
      double x = coord.first; double y = coord.second;
      return 2*(x*x + y*y);
    }
    
    double boundary_fn( const std::pair<double,double>& coord ) {
      double x = coord.first; double y = coord.second;
      return x*x*y*y;
    }
    
  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main(int argc, char *argv[])
{
  PetscSession::getInstance(argc,argv);

  // Number of points in a square mesh.
  unsigned Nx = 5;
  unsigned Ny = 5;
  
  DenseVector<double> X = Utility::uniform_node_vector( -2.0, 2.0, Nx );
  DenseVector<double> Y = Utility::uniform_node_vector( -1.0, 1.0, Ny );
  TwoD_Node_Mesh<double> mesh(X,Y,2);
  double dx2 = pow(X[1]-X[0],2);
  double dy2 = pow(Y[1]-Y[0],2);
  
  cout << "\n";
  cout << "=== Poisson: Cartesian geometry =====================\n";
  cout << "    Solving in [-2 , 2] x [-1 , 1]      \n";
  cout << "    Using a " << Nx << "x" << Ny << " mesh\n";
  cout << "\n";

  SparseMatrix<double> A(Nx*Ny,Nx*Ny);
  DenseVector<double> B(Nx*Ny,0.0);
  {
    // LHS
    unsigned i(0);
    for (unsigned j = 0; j < Ny; ++j) {
      A(j,j) = 1.0;
      B[j] = Example::boundary_fn(mesh.coord(i,j));
    }
  }
  // step through x nodes
  for (unsigned i = 1; i < Nx-1; ++i) {
    // bottom
    unsigned j = 0;
    A(i*Ny+j,i*Ny+j) = 1.0;
    B[i*Ny+j] = Example::boundary_fn(mesh.coord(i,j));
    for ( j = 1; j < Ny-1; ++j) {
      // interior nodes
      A(i*Ny+j,i*Ny+j-1) = 1./dy2;
      A(i*Ny+j,i*Ny+j) = -2./dx2-2./dy2;
      A(i*Ny+j,i*Ny+j+1) = 1./dy2;
      A(i*Ny+j,i*Ny+j+Ny) = 1./dx2;
      A(i*Ny+j,i*Ny+j-Ny) = 1./dx2;
      B[i*Ny+j] = Example::source_fn(mesh.coord(i,j));
    }
    // top
    j = Ny-1;
    A(i*Ny+j,i*Ny+j) = 1.0;
    B[i*Ny+j] = Example::boundary_fn(mesh.coord(i,j));
  }
  {
    // RHS
    unsigned i(Nx-1);
    for (unsigned j = 0; j < Ny; ++j) {
      A(i*Ny+j,i*Ny+j) = 1.0;
      B[i*Ny+j] = Example::boundary_fn(mesh.coord(i,j));
    }
  }

  //A.dump();
  //B.dump();

  SparseLinearSystem<double> system(&A,&B,"petsc");
  system.solve();

  for (unsigned i = 0; i < Nx; ++i) {
    for (unsigned j = 0; j < Ny; ++j) {
      const double x = mesh.coord(i,j).first; 
      const double y = mesh.coord(i,j).second;
      // solution is quadratic in x,y so should be accurate on
      // any mesh
      mesh(i,j,0)=B[i*Nx+j];
      mesh(i,j,1)=mesh(i,j,0) - pow(x*y,2.0);
    }
  }

  mesh.dump_gnu("./DATA/mesh.dat");
  
  bool failed(true);
  const double tol(1.e-14);
  double abserror = mesh.max(1);
  cout << "error = " << abserror << "\n";

  if ( abserror < tol ) failed = false;

  if ( failed )
  {
    cout << "Poisson solver failed to give exact solution.\n";
    cout << "error = " << abserror << "\n";
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

}
