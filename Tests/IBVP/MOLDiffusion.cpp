/// Quick Method Of Lines test.

#include <IVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    unsigned N(201);
    
    /// Define the equations by inheriting Equation base class
    class Diff_equation : public Equation<double>
    {
    public:

      /// Constructor
      Diff_equation() : Equation<double>( N ) {};

      /// We implement the equation as N first-order ODEs
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        // BC at y=0
        f[0] = 0;
        for (unsigned i = 1; i < N-1; ++i) {
          // heat equation
          f[i] = ( z[i-1]-2*z[i]+z[i+1] )/(dy*dy);
        }
        // BC at y=1
        f[N-1] = 0;
      }

      /// The usual 3 parameters of the Lorenz eqns
      double dy;

    };
  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== MOL: RK integration of the Diffusion equation ===\n";
  cout << "\n";

  DenseVector<double> nodes = Utility::uniform_node_vector(0.0,1.0,Example::N);
  DenseVector<double> u( Example::N, 0.0 );
  // initialise
  u[ 0 ] = 0.0;
  u[ Example::N-1 ] = 0.0;
  for (unsigned i = 1; i < Example::N-1; ++i) {
    double y = nodes[i];
    u[ i ] = y*(1-y);
  }
  
  const int max_num_of_steps( 100000 );
  // set up the problem.
  Example::Diff_equation problem;
  problem.dy = nodes[1]-nodes[0];
  //cout << "dy^2=" << problem.dy*problem.dy << " dt = " << 1.0/max_num_of_steps << " dy^2/dt = " << problem.dy*problem.dy*max_num_of_steps << "\n";
  
  // Construct an ODE from the problem.
  ODE_IVP<double> ode( &problem, 0.0, 1.0, max_num_of_steps );
  ode.store_every() = 100;
  u = ode.shoot4( u );
  OneD_Node_Mesh<double> soln = ode.get_mesh();
  
  // output and check
  const double tol(1.e-4);
  std::size_t actual_num_of_steps = soln.get_nnodes();
  double max_diff(0.0);
  for (std::size_t i=0; i<actual_num_of_steps; ++i) {
    double time = soln.coord(i);
    // evaluate analytical Fourier series solution at y = 0.5
    double y = 0.5;
    double u( 0.0 );
    int en( -1 );
    double correction( 0.0 );      
    do
      {
        en += 2;
        correction = 8 / ( std::pow( en * M_PI, 3 ) )
          * std::exp( -std::pow( en * M_PI, 2 ) * time ) * std::sin( en * M_PI * y );
        u += correction;
      }
    while ( std::abs( correction ) > tol / 10. );
    // examine the difference between the numerical and series solutions
    max_diff = std::abs( u - soln(i,(Example::N-1)/2));
  }

  bool failed = true;
  if ( abs( max_diff ) < tol )
  {
    failed = false;
  }

  if ( failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }  
  
}
