/// \file HYP_2D_acoustic_NV.cpp
/// \ingroup Examples
/// \ingroup HYP_2D
/// A linear acoustic pulse propagating towards a cylindrical
/// region of higher (four times) bulk modulus.

#include <TwoD_HYP_bundle.h>

enum { p, u, v };

namespace CppNoddy
{
  namespace Example
  {
    /// bulk modulus
    double K_small( 0.1 );
    double K_big( 1.0 );
    /// lower density
    double rho_small( 1.0 );

    /// A smoothed 1-D top hat function for the change in material properties
    double hat( const double &x, const double &xr, const double &xl, const double &sharpness )
    {
      return ( std::tanh( sharpness * ( x - xl ) ) - std::tanh( sharpness * ( x - xr ) ) ) / 2.0;
    }

    double rho( const double &x, const double &y )
    {
      return rho_small;
    }

    double K( const double &x, const double &y )
    {
      return K_small + hat( x, 0.2, -0.2, 50. ) * hat( y, 0.2, -0.2, 50. ) * ( K_big - K_small );
    }

    /// Define the system
    class Acoustic_2d : public TwoD_Hyperbolic_System
    {

    public:

      // A 2D advection equation
      Acoustic_2d() : TwoD_Hyperbolic_System( 3 )
      {}

      /// Define the vector flux in the x direction
      void flux_fn_x( const DenseVector<double>& x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ p ] = K( x[0], x[1] ) * q[ u ];
        f[ u ] = q[ p ] / rho( x[0], x[1] );
        f[ v ] = 0.0;
      }

      /// Define the vector flux in the y direction
      void flux_fn_y( const DenseVector<double>& x, const DenseVector<double> &q, DenseVector<double> &g ) const
      {
        g[ p ] = K( x[0], x[1] ) * q[ v ];
        g[ u ] = 0.0;
        g[ v ] = q[ p ] / rho( x[0], x[1] );
      }

      /// Bound the wave speed
      void max_charac_speed( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &c ) const
      {
        // maximum wave speed
        c[ 0 ] = c[ 1 ] = sqrt( K_big / rho_small );
      }

      void source_fn( const DenseVector<double>& x, const DenseVector<double>& q, DenseVector<double>& r ) const
      {
        const double delta( 1.e-6 );
        r[ p ] = + q[ u ] * ( K( x[0] + delta, x[1] ) - K( x[0], x[1] ) ) / delta
                 + q[ v ] * ( K( x[0], x[1] + delta ) - K( x[0], x[1] ) ) / delta;
        r[ u ] = - q[ p ] * ( ( rho( x[0] + delta, x[1] ) - rho( x[0], x[1] ) ) / delta ) / std::pow( rho( x[0], x[1] ), 2 );
        r[ v ] = - q[ p ] * ( ( rho( x[0], x[1] + delta ) - rho( x[0], x[1] ) ) / delta ) / std::pow( rho( x[0], x[1] ), 2 );
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, const double &y, DenseVector<double> &q )
    {
      if ( ( x > -0.4 ) && ( x < -0.35 ) )
      {
        q[ p ] = 1.0;
        q[ u ] = q[ p ] / sqrt( rho( x, y ) * K( x, y ) );
        q[ v ] = 0.0;
      }
      else
      {
        q[ p ] = q[ u ] = q[ v ] = 0.0;
      }
    }
  } //end Example namespace

} //end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  std::string filename_stub( "./DATA/HYP_2D_acoustic" );
  // define the domain/mesh
  const double west =  -0.5;
  const double east = 0.5;
  const double south = -0.5;
  const double north = 0.5;
  const unsigned Nx = 201;
  const unsigned Ny = 201;
  DenseVector<double> faces_x = Utility::power_node_vector( west, east, Nx, 1.0 );
  DenseVector<double> faces_y = Utility::power_node_vector( south, north, Ny, 1.0 );

  Example::Acoustic_2d conservative_problem;
  TwoD_TVDLF_Mesh Acoustic_mesh( faces_x, faces_y, &conservative_problem, Example::Q_init );
  Acoustic_mesh.set_limiter( 1 );
  Acoustic_mesh.dump_nodes_x( "./DATA/HYP_2D_acoustic.xnodes" );
  Acoustic_mesh.dump_nodes_y( "./DATA/HYP_2D_acoustic.ynodes" );
  std::ofstream kvalues;
  std::string filename( "./DATA/HYP_2D_impedance.dat" );
  kvalues.open( filename.c_str() );
  for ( unsigned i = 1; i < Nx; ++i )
  {
    for ( unsigned j = 1; j < Ny; ++j )
    {
      double x( ( faces_x[ i - 1 ] + faces_x[ i ] ) / 2.0 );
      double y( ( faces_y[ j - 1 ] + faces_y[ j ] ) / 2.0 );
      kvalues << Example::K( x, y ) << "\n";
      kvalues.flush();
    }
  }

  int file_counter( 1 );
  int loop_counter( 0 );
  const double t_end = 10.;
  do
  {
    Acoustic_mesh.update( 0.49, std::abs( Acoustic_mesh.get_time() - t_end ) );
    if ( loop_counter % 4 == 0 )
    {
      cout << " *** current time = " << Acoustic_mesh.get_time() << "\n";
      Acoustic_mesh.dump_data( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      Acoustic_mesh.dump_gnu( filename_stub + Utility::stringify( file_counter ) + "_gnu.dat" );
      file_counter += 1;
    }
    ++loop_counter;
  }
  while ( Acoustic_mesh.get_time() < t_end );

}
