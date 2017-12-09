/// \file HYPSodsShocktube.cpp
/// \ingroup Tests
/// \ingroup HYP_1D
/// Solving the 1D Euler problem for gas dynamics
/// \f[ \rho_t + \left ( m \right )_x = 0 \f]
/// \f[ m_t + \left ( \rho u^2 + p \right )_x = 0 \f]
/// \f[ E_t + \left ( u ( E + p ) \right )_x = 0 \f]
/// where
/// \f[ u = m / \rho \f]
/// \f[ p = ( \gamma - 1 ) ( E - \frac12 \rho u^2 ) \f]
/// and \f$ \gamma = 1.4 \f$.
/// The initial conditions correspond to Sod's problem
/// \f[ (\rho, m, E ) = (  1, 0, 2.5 ) \f] for \f$ x < \frac12 \f$
/// \f[ (\rho, m, E ) = ( 0.125, 0, 0.25 ) \f] for \f$ x > \frac12 \f$.

#include <OneD_HYP_bundle.h>

enum { rho, m, E };

namespace CppNoddy
{
  namespace Example
  {
    /// Adiabatic index -- ratio of specific heats
    double gamma( 1.4 );

    /// Define the system
    class Euler_1d : public OneD_Hyperbolic_System
    {

    public:

      /// One dimemsional Euler problem, so 3rd order
      Euler_1d() : OneD_Hyperbolic_System( 3 )
      {}

      /// Define the vector flux
      void flux_fn( const double &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        double u = q[ m ] / q[ rho ];
        double p = ( gamma - 1. ) * ( q[ E ] - 0.5 * q[ rho ] * u * u );
        f[ rho ] = q[ m ];
        f[ m ]   = q[ rho ] * u * u + p;
        f[ E ]   = u * ( q[ E ] + p );
      }

      /// Bound the shock speed
      double max_charac_speed( const DenseVector<double> &q ) const
      {
        // sound speeds
        double u = q[ m ] / q[ rho ];
        double p = ( gamma - 1. ) * ( q[ E ] - 0.5 * q[ rho ] * u * u );
        double c = sqrt( gamma *  p / q[ rho ] );
        // maximum shock speed
        return std::max(  u + c, u - c );
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      if ( x < 0.5 )
      {
        q[ rho ] = 1.0;
        q[ m ] = 0.0;
        q[ E ] = 2.5; // => P = 1.0
      }
      else
      {
        q[ rho ] = 0.125;
        q[ m ] = 0.0;
        q[ E ] = 0.25; // => P = 0.1
      }
    }
  } //end Example namespace
} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== Hyperbolic: 1D Euler gasdynamics shocktube  =====\n";
  cout << "\n";

  // define the domain/mesh
  const double left =  0.0;
  const double right = 1.0;
  const unsigned N = 400;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  double t = 0.0;

  Example::Euler_1d conservative_problem;
  OneD_TVDLF_Mesh Euler_mesh( faces_x, &conservative_problem, Example::Q_init );
  Euler_mesh.set_limiter( 0 );

  double I1 = Euler_mesh.integrate()[0];
  int loop_counter( 0 );
  int file_counter( 0 );

  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_shocktube_sod";
  TrackerFile my_file( 10 );
  OneD_Node_Mesh<double> soln = Euler_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    if ( loop_counter % 10 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Euler_mesh.get_soln( );
      my_file.update();
      file_counter += 1;
    }
    t += Euler_mesh.update( 0.499 );
    ++loop_counter;
  }
  while ( ( t < 0.2 ) && ( loop_counter < 1000 ) );

  double I2 = Euler_mesh.integrate()[0];
  // These are clawpack results (to 4dp) with N=400 and minmod for the
  // density of the two density steps ...
  soln = Euler_mesh.get_soln( );
  double E1 = std::abs( soln.get_interpolated_vars( 0.575 )[0] - 0.4264 );
  double E2 = std::abs( soln.get_interpolated_vars( 0.75 )[0] - 0.2656 );
  // check against the clawpack data, and that the mass is conserved.
  if ( ( E1 > 1.e-3 ) || ( E2 > 1.e-3 ) || ( std::abs( I1 - I2 ) > 1.e-8 ) || ( loop_counter >= 1000 ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << E1 << " " << E2 << " " << std::abs( I1 - I2 ) << " " << loop_counter << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

} // end of main()
