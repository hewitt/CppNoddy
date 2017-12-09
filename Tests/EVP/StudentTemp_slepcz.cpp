/// TEMP FILE

#include <EVP_bundle.h>
#include <BVP_bundle.h>
#include <SparseLinearSystem.h>

using namespace CppNoddy;
using namespace std;

enum {q,phi};

int main()
{
  cout << "\n";
  cout << "=== EVP: Temporal spectra of the Orr-Sommerfeld eqn =\n";
  cout << "===  with a matrix problem assembled by hand.\n";
  cout << "\n";


  // discretise with these many nodal points
  const std::size_t nnodes( 201 );
  // domain boundaries
  const double left = -1.0;
  const double right = 1.0;
  DenseVector<double> y_nodes = Utility::uniform_node_vector(left, right, nnodes);
  const double d = y_nodes[1]-y_nodes[0];
  // GUESS for the local refinement
  OneD_Node_Mesh<D_complex> guess( y_nodes, 2 );
  D_complex C(0.0,0.0);

  // we'll solve as TWO second order problems
  const std::size_t N( 2 * nnodes );

  // streamwise wavenumber and Reynolds number
  double alpha ( 1.02 );
  double Re ( 5772.2 );
  const D_complex I( 0.0, 1.0 );

  {
    // matrices for the EVP, initialised with zeroes
    DenseMatrix<D_complex> a( N, N, 0.0 );
    DenseMatrix<D_complex> b( N, N, 0.0 );
    // boundary conditions at the left boundary
    a( 0, 0 ) = 1.0;           // phi( left ) = 0
    a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
    a( 1, 2 ) = 2.0 / d;
    a( 1, 4 ) = -0.5 / d;
    // fill the interior nodes
    for ( std::size_t i = 1; i <= nnodes - 2; ++i )
    {
      // position in the channel
      const double y = y_nodes[i];
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
    catch ( std::runtime_error )
    {
      cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
      return 1;
    }
    // tag any eigenvalues with imaginary part > -0.01 -- there is only 1
    system.set_shift( D_complex( 0.0, -0.01 ) );
    system.tag_eigenvalues_upper( + 1 );
    lambdas = system.get_tagged_eigenvalues();
    lambdas.dump();
    DenseMatrix<D_complex> evecs( system.get_tagged_eigenvectors() );
    // store the poor mesh result as a guess
    for ( unsigned i=0; i<nnodes; ++i )
    {
      guess(i,0)=evecs(0,2*i+0);
      guess(i,1)=evecs(0,2*i+1);
    }
    guess.scale( 1.0/guess(nnodes-1,1) );
    C = lambdas[0];
  }
  guess.dump_gnu("vec.dat");

  PetscInitialize(NULL,NULL,(char*)0,(char*)0);
  double max_residual(1.0);
  int iteration(1);
  Re = 5800;//772.2;

  /*
  do
  {
    // matrices for the EVP, initialised with zeroes
    SparseMatrix<D_complex> a( N+1, N+1 );
    // RHS residual vector
    DenseVector<D_complex> B( N+1, 0.0 );

    // boundary conditions at the left boundary
    a( 0, 0 ) = 1.0;           // phi( left ) = 0
    B[0] = - (guess(0,0));
    a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
    a( 1, 2 ) = 2.0 / d;
    a( 1, 4 ) = -0.5 / d;
    B[1] = -( -1.5*guess(0,0)+2.0*guess(0,0)-0.5*guess(0,0) )/d;
    // fill the interior nodes
    for ( std::size_t i = 1; i <= nnodes - 2; ++i )
    {
      // position in the channel
      const double y = y_nodes[i];
      // Poiseuille flow profile
      const double U = ( 1.0 - y * y );
      const double Udd = -2.0;

      // the first quation at the i'th nodal point
      std::size_t row = 2 * i;
      a( row, row ) = -2.0 / ( d * d ) - alpha * alpha;
      a( row, row - 2 ) = 1.0 / ( d * d );
      a( row, row + 2 ) = 1.0 / ( d * d );
      a( row, row + 1 ) = -1.0;
      B[row] = -( guess(i,1) - (guess(i-1,0)-2.0*guess(i,0)+guess(i+1,0))/(d*d) + alpha*alpha*guess(i,0));

      row += 1;
      // the second equation at the i'th nodal point
      a( row, row ) = -2.0 / ( d * d ) - alpha * alpha - I * alpha * Re * U + I*alpha*Re*C;
      a( row, row - 2 ) = 1.0 / ( d * d );
      a( row, row + 2 ) = 1.0 / ( d * d );
      a( row, row - 1 ) = I * alpha * Re * Udd;
      a( row, N ) = I*alpha*Re*guess(i,1);
      B[row] = -( (guess(i-1,1)-2.0*guess(i,1)+guess(i+1,1))/(d*d) - alpha*alpha*guess(i,1) - I*alpha*Re*U*guess(i,1)
          + I*alpha*Re*Udd*guess(i,0) + I*alpha*Re*C*guess(i,1) );
    }
    // boundary conditions at right boundary
    a( N - 2, N - 2 ) = 1.5 / d;
    a( N - 2, N - 4 ) = -2.0 / d;
    a( N - 2, N - 6 ) = 0.5 / d; // psi'( right ) = 0
    B[N-2] = -( 1.5*guess(nnodes-1,0)-2.0*guess(nnodes-2,0)+0.5*guess(nnodes-3,0) )/d;
    a( N - 1, N - 2 ) = 1.0;     // psi( right ) = 0
    B[N-1] = -( guess(nnodes-1,0) );
    a( N, N-1 ) = 1.0;
    B[N] = 1.0 - guess(nnodes-1,1);

    //a.dump();
    //B.dump();

    // a vector for the eigenvalues
    max_residual = B.inf_norm();
    //cout << " inf_norm(B) = " << B.inf_norm() << "\n";
    SparseLinearSystem<D_complex> linsys( &a, &B, "petsc" );
    linsys.solve();
    // a vector for the eigenvalues
    cout << "iteration = " << iteration << " inf_norm(correction) = " << B.inf_norm() << " current c = " << C << "\n";
    iteration += 1;


    for ( std::size_t i = 0; i < nnodes; ++i )
    {
      guess(i,0) += B[i*2+0];
      guess(i,1) += B[i*2+1];
    }
    C += B[N];
  } while (max_residual > 1.e-8);
  */

top:

  //C += D_complex(.01,.01);
  do
  {
    D_complex oldf(0.0);
    D_complex newf(0.0);
    {
      // matrices for the EVP, initialised with zeroes
      SparseMatrix<D_complex> a( N, N );
      // RHS residual vector
      DenseVector<D_complex> B( N, 0.0 );

      // boundary conditions at the left boundary
      a( 0, 0 ) = 1.0;           // phi( left ) = 0
      B[0] = 1.0;
      a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
      a( 1, 2 ) = 2.0 / d;
      a( 1, 4 ) = -0.5 / d;
      // fill the interior nodes
      for ( std::size_t i = 1; i <= nnodes - 2; ++i )
      {
        // position in the channel
        const double y = y_nodes[i];
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
        a( row, row ) = -2.0 / ( d * d ) - alpha*alpha - I*alpha*Re*U + I*alpha*Re*C;
        a( row, row - 2 ) = 1.0 / ( d * d );
        a( row, row + 2 ) = 1.0 / ( d * d );
        a( row, row - 1 ) = I * alpha * Re * Udd;
      }
      // boundary conditions at right boundary
      a( N - 2, N - 2 ) = 1.5 / d;
      a( N - 2, N - 4 ) = -2.0 / d;
      a( N - 2, N - 6 ) = 0.5 / d; // psi'( right ) = 0
      a( N - 1, N - 2 ) = 1.0;     // psi( right ) = 0

      // a vector for the eigenvalues
      max_residual = B.inf_norm();
      //cout << " inf_norm(B) = " << B.inf_norm() << "\n";
      SparseLinearSystem<D_complex> linsys( &a, &B, "petsc" );
      linsys.solve();

      D_complex sum(0.0);
      for ( std::size_t i = 1; i < nnodes; ++i )
      {
        sum += 0.5*(B[(i-1)*2]*B[(i-1)*2]+B[i*2]*B[i*2])*(y_nodes[i]-y_nodes[i-1]);
        //sum += B[i*2+0]*(y_nodes[i]-y_nodes[i-1]);
      }
      // a vector for the eigenvalues
      //cout << "C = " << C << " inv_integral = " << 1./sum << "\n";
      oldf = sqrt(1./sum);
    }


    double delta( 1.e-6 );
    C += delta;


    {
      // matrices for the EVP, initialised with zeroes
      SparseMatrix<D_complex> a( N, N );
      // RHS residual vector
      DenseVector<D_complex> B( N, 0.0 );

      // boundary conditions at the left boundary
      a( 0, 0 ) = 1.0;           // phi( left ) = 0
      B[0] = 1.0;
      a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
      a( 1, 2 ) = 2.0 / d;
      a( 1, 4 ) = -0.5 / d;
      // fill the interior nodes
      for ( std::size_t i = 1; i <= nnodes - 2; ++i )
      {
        // position in the channel
        const double y = y_nodes[i];
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
        a( row, row ) = -2.0 / ( d * d ) - alpha*alpha - I*alpha*Re*U + I*alpha*Re*C;
        a( row, row - 2 ) = 1.0 / ( d * d );
        a( row, row + 2 ) = 1.0 / ( d * d );
        a( row, row - 1 ) = I * alpha * Re * Udd;
      }
      // boundary conditions at right boundary
      a( N - 2, N - 2 ) = 1.5 / d;
      a( N - 2, N - 4 ) = -2.0 / d;
      a( N - 2, N - 6 ) = 0.5 / d; // psi'( right ) = 0
      a( N - 1, N - 2 ) = 1.0;     // psi( right ) = 0

      // a vector for the eigenvalues
      max_residual = B.inf_norm();
      //cout << " inf_norm(B) = " << B.inf_norm() << "\n";
      SparseLinearSystem<D_complex> linsys( &a, &B, "petsc" );
      linsys.solve();

      D_complex sum(0.0);
      for ( std::size_t i = 1; i < nnodes; ++i )
      {
       sum += 0.5*(B[(i-1)*2]*B[(i-1)*2]+B[i*2]*B[i*2])*(y_nodes[i]-y_nodes[i-1]);
       //sum += B[i*2+0]*(y_nodes[i]-y_nodes[i-1]);
      }
      // a vector for the eigenvalues
      //cout << "C = " << C << " inv_integral = " << 1./sum << "\n";
      newf = sqrt(1./sum);
    }

    //cout << " oldf = " << oldf << " newf = " << newf << " f' = " << (newf-oldf)/1.e-4 << "\n";
    //cout << " delta = " << -oldf/((newf-oldf)/1.e-8) << "\n";

    max_residual = abs(0.5*(oldf+newf));

    cout << "residual = " << max_residual << "\n";
    C += -0.5*(oldf+newf)/((newf-oldf)/delta);
    //cout << "C = " << C << "\n";

    //return 0;
  } while ( max_residual > 1.e-8 );

  cout << Re << " " << C.real() << " " << C.imag() << "\n";
  Re += 10.0;
  if ( Re < 20000 )
  {
    goto top;
  }
}
