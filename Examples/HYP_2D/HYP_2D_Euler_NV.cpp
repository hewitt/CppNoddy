/// \file HYP_2D_Euler_NV.cpp
/// \ingroup Examples
/// \ingroup HYP_2D
/// Solving the Euler equations for a compressible gas
/// with Sod-like initial conditions in a
/// rectangular box \f$ [0,0.3] \times [0,0.3]\f$,
/// with no-momentum flux through the walls.
/// The initial state is one of no motion but a step
/// change in density and pressure along the line
/// \f$ x + y = 0.15 \f$.

#include <TwoD_HYP_bundle.h>

enum { rho, mx, my, E };

namespace CppNoddy
{
  namespace Example
  {
    /// Adiabatic index -- ratio of specific heats
    double gamma( 1.4 );

    /// Define the system
    class Euler_2d : public TwoD_Hyperbolic_System
    {

    public:

      /// One dimemsional Euler problem, so 3rd order
      Euler_2d() : TwoD_Hyperbolic_System( 4 )
      {}

      /// Define the vector flux
      void flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        double u = q[ mx ] / q[ rho ];
        double v = q[ my ] / q[ rho ];
        double p = ( gamma - 1. ) * ( q[ E ] - 0.5 * q[ rho ] * ( u * u + v * v ) );
        //
        f[ rho ] = q[ mx ];
        f[ mx ]  = q[ rho ] * u * u + p;
        f[ my ]  = q[ rho ] * u * v;
        f[ E ]   = u * ( q[ E ] + p );
      }
      void flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &g ) const
      {
        double u = q[ mx ] / q[ rho ];
        double v = q[ my ] / q[ rho ];
        double p = ( gamma - 1. ) * ( q[ E ] - 0.5 * q[ rho ] * ( u * u + v * v ) );
        //
        g[ rho ] = q[ my ];
        g[ mx ]  = q[ rho ] * u * v;
        g[ my ]  = q[ rho ] * v * v + p;
        g[ E ]   = v * ( q[ E ] + p );
      }

      /// Bound the wave speed
      void max_charac_speed( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &c ) const
      {
        // maximum wave speed
        const double u = q[ mx ] / q[ rho ];
        const double v = q[ my ] / q[ rho ];
        const double px = ( gamma - 1. ) * ( q[ E ] - 0.5 * q[ rho ] * u * u );
        const double cx = sqrt( gamma *  px / q[ rho ] );
        const double py = ( gamma - 1. ) * ( q[ E ] - 0.5 * q[ rho ] * v * v );
        const double cy = sqrt( gamma *  py / q[ rho ] );
        // maximum shock speed
        c[ 0 ] = std::max( u + cx, u - cx );
        c[ 1 ] = std::max( v + cy, v - cy );
      }

      /// edge conditions
      std::vector<bool> edge_values( const int& face_index, const DenseVector<double>& x, DenseVector<double>& q ) const
      {
        std::vector<bool> comps( q.size(), false );
        switch ( face_index )
        {
        case 0:
          q[ my ] = 0.0;
          comps[ my ] = true;
          break;
        case 1:
          q[ mx ] = 0.0;
          comps[ mx ] = true;
          break;
        case 2:
          q[ my ] = 0.0;
          comps[ my ] = true;
          break;
        case 3:
          q[ mx ] = 0.0;
          comps[ mx ] = true;
          break;
        }
        return comps;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, const double &y, DenseVector<double> &q )
    {
      const double xi( x + y );
      const double factor( 500.0 );
      const double rho_1( 1.0 );
      const double rho_2( 0.125 );
      const double E_1( 1.0 / ( gamma - 1.0 ) );
      const double E_2( 0.14 / ( gamma - 1.0 ) );
      q[ rho ] = ( rho_1 + rho_2 ) * 0.5 + ( rho_1 - rho_2 ) * 0.5 * std::tanh( factor * ( xi - 0.15 ) );
      q[ E ] = ( E_1 + E_2 ) * 0.5 + ( E_1 - E_2 ) * 0.5 * std::tanh( factor * ( xi - 0.15 ) );
    }
  } //end Example namespace
} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  std::string filename_stub( "./DATA/HYP_2D_Euler" );
  // define the domain/mesh
  const double west = 0.0;
  const double east = 0.3;
  const double south = 0.0;
  const double north = 0.3;
  const unsigned Nx = 201;
  const unsigned Ny = 201;
  // non-uniform mesh, with more points at the focussing corner
  DenseVector<double> faces_x = Utility::power_node_vector( west, east, Nx, 1.1 );
  DenseVector<double> faces_y = Utility::power_node_vector( south, north, Ny, 1.1 );

  Example::Euler_2d conservative_problem;
  TwoD_TVDLF_Mesh Euler_2d_mesh( faces_x, faces_y, &conservative_problem, Example::Q_init );
  Euler_2d_mesh.set_limiter( 1 );
  Euler_2d_mesh.dump_nodes_x( "./DATA/HYP_2D_Euler.xnodes" );
  Euler_2d_mesh.dump_nodes_y( "./DATA/HYP_2D_Euler.ynodes" );

  int file_counter( 1 );
  int loop_counter( 0 );
  const double t_end = 1.;
  do
  {
    if ( loop_counter % 10 == 0 )
    {
      cout << " *** " << file_counter << " current time = " << Euler_2d_mesh.get_time() << "\n";
      Euler_2d_mesh.dump_data( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      file_counter += 1;
    }
    Euler_2d_mesh.update( 0.49 );//, std::abs( Euler_2d_mesh.get_time() - t_end ) );
    loop_counter += 1;
  }
  while ( Euler_2d_mesh.get_time() < t_end );


} // end of main()
