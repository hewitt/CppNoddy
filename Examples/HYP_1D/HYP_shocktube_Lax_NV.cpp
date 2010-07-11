/// \file HYP_shocktube_Lax_NV.cpp
/// \ingroup Examples
/// \ingroup HYP_1D
/// Solving the 1D Euler problem for gas dynamics
/// \f[ \rho_t + \left ( m \right )_x = 0 \f]
/// \f[ m_t + \left ( \rho u^2 + p \right )_x = 0 \f]
/// \f[ E_t + \left ( u ( E + p ) \right )_x = 0 \f]
/// where
/// \f[ u = m / \rho \f]
/// \f[ p = ( \gamma - 1 ) ( E - \frac12 \rho u^2 ) \f]
/// and \f$ \gamma = 1.4 \f$.
/// The initial conditions correspond to Lax's problem
/// \f[ (\rho, m, E ) = ( 0.445, 0.311, 8.928 ) \f] for \f$ x < \frac12 \f$
/// \f[ (\rho, m, E ) = ( 0.5, 0, 1.428 ) \f] for \f$ x > \frac12 \f$.

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

      /// edge conditions
      std::vector<bool> edge_values( const int face_index, const double& x, DenseVector<double>& q ) const
      {
        // inflow at all boundaries?
        std::vector<bool> inflow( q.size(), true );
        // x doesn't matter since the conditions are fixed
        if ( face_index == -1 )
        {
          q[ rho ] = 0.445;
          q[ m ] = 0.445 * 0.698;
          q[ E ] = 3.528 / ( gamma - 1. ) + 0.5 * 0.445 * 0.698 * 0.698; // => P = 3.528
        }
        else
        {
          q[ rho ] = 0.5;
          q[ m ] = 0.0;
          q[ E ] = 0.571 / ( gamma - 1. ); // => P = 0.571
        }
        return inflow;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      if ( x < 0.5 )
      {
        q[ rho ] = 0.445;
        q[ m ] =  0.445 * 0.698;
        q[ E ] = 3.528 / ( gamma - 1. ) + 0.5 * 0.445 * 0.698 * 0.698;
      }
      else
      {
        q[ rho ] = 0.5;
        q[ m ] = 0.0;
        q[ E ] = 0.571 / ( gamma - 1. );
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
  cout << "    Lax's initial conditions \n";
  cout << "\n";

  // define the domain/mesh
  const double left =  0.0;
  const double right = 1.0;
  const unsigned N = 401;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  double t = 0.0;

  Example::Euler_1d conservative_problem;
  OneD_TVDLF_Mesh Euler_mesh( faces_x, &conservative_problem, Example::Q_init );
  Euler_mesh.set_limiter( 0 );

  int loop_counter( 2 );
  int file_counter( 1 );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_shocktube_Lax";
  TrackerFile my_file( 10 );
  OneD_Node_Mesh<double> soln = Euler_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    t += Euler_mesh.update( 0.49 );
    std::cout << t << "\n";
    if ( loop_counter % 10 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Euler_mesh.get_soln( );
      my_file.update();
      file_counter += 1;
    }
    ++loop_counter;
  }
  while ( ( t < 0.1 ) && ( loop_counter < 1000 ) );

  my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
  my_file.update();


} // end of main()
