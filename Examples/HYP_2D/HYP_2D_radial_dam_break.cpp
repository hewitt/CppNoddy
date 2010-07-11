/// \file HYP_2D_radial_dam_break.cpp
/// \ingroup Examples
/// \ingroup HYP_2D
/// A radial dam-break problem for the 2D shallow water equations.
/// The flow height at a single point is compared with the result
/// obtained from CLAWPACK.

#include <string>

#include <Utility.h>
#include <TwoD_HYP_bundle.h>

enum { h, hu, hv };

namespace CppNoddy
{
  namespace Example
  {
    /// gravitational acceleration
    double g( 1.0 );
    /// initial hump amplitude
    double A( 1.0 );
    /// Set the initial state of the system
    void Q_init( const double &x, const double &y, DenseVector<double> &q )
    {
      double r( std::sqrt( x*x + y*y ) - 0.5 );
      q[ h ] = 1 + 0.5 * ( 1 - std::tanh( 500 * r ) );
      q[ hu ] = 0;
    }

    /// Define the system
    class Shallow_2d_rad : public TwoD_Hyperbolic_System
    {

    public:

      /// One dimemsional constant coefft acoustic problem
      Shallow_2d_rad() : TwoD_Hyperbolic_System( 3 )
      {}

      /// Define the vector flux
      void flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ h ] = q[ hu ];
        f[ hu ] = q[ hu ] * q[ hu ] / q[ h ] + g * q[ h ] * q[ h ] / 2.;
        f[ hv ] = q[ hu ] * q[ hv ] / q[ h ];
      }

      void flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ h ] = q[ hv ];
        f[ hu ] = q[ hu ] * q[ hv ] / q[ h ];
        f[ hv ] = q[ hv ] * q[ hv ] / q[ h ] + g * q[ h ] * q[ h ] / 2.;
      }

      // just for fun, we'll use analytic Jacobians
      void Jac_flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseMatrix<double> &J ) const
      {
        J( 0, hu ) = 1;
        J( 1, h  ) = - q[ hu ] * q[ hu ] / ( q[ h ] * q[ h ] ) + g * q[ h ];
        J( 1, hu ) = 2 * q[ hu ] / q[ h ];
        J( 2, h  ) = - q[ hu ] * q[ hv ] / ( q[ h ] * q[ h ] );
        J( 2, hu ) = q[ hv ] / q[ h ];
        J( 2, hv ) = q[ hu ] / q[ h ];
      }

      void Jac_flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseMatrix<double> &J ) const
      {
        J( 0, hv ) = 1;
        J( 1, h  ) = - q[ hu ] * q[ hv ] / ( q[ h ] * q[ h ] );
        J( 1, hu ) = q[ hv ] / q[ h ];
        J( 1, hv ) = q[ hu ] / q[ h ];
        J( 2, h  ) = - q[ hv ] * q[ hv ] / ( q[ h ] * q[ h ] ) + g * q[ h ];
        J( 2, hv ) = 2 * q[ hv ] / q[ h ];
      }

      /// Bound the wave speed
      void max_charac_speed( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &c ) const
      {
        // wave speed
        double cc = sqrt( g * q[ h ] );
        // flow speed
        double U = q[ hu ] / q[ h ];
        double V = q[ hv ] / q[ h ];
        // maximum shock speed
        c[ 0 ] = std::max( std::abs( U + cc ), std::abs( U - cc ) );
        c[ 1 ] = std::max( std::abs( V + cc ), std::abs( V - cc ) );
      }
      
    };

  } //end Example namespace

} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  std::string filename_stub( "./DATA/HYP_2D_rad_dam" );

  cout << "\n";
  cout << "=== Hyperbolic: 2D radial dam-break in Cartesians ===\n";
  cout << "\n";

  // define the domain/mesh
  const double west = -1;
  const double east = 1;
  const double south = -1;
  const double north = 1;
  const unsigned N = 151;
  DenseVector<double> faces_x = Utility::power_node_vector( west, east, N, 1.0 );
  DenseVector<double> faces_y = Utility::power_node_vector( south, north, N, 1.0 );

  Example::Shallow_2d_rad conservative_problem;
  TwoD_TVDLF_Mesh Shallow_2d_mesh( faces_x, faces_y, &conservative_problem, Example::Q_init );
  Shallow_2d_mesh.set_limiter( 0 );

  int file_counter( 1 );
  const double t_end = 0.5;
  DenseVector<double> x1(2, 0.0); x1[0]=0.4; x1[1]=0.1;
  DenseVector<double> x2(2, 0.0); x2[0]=0.1; x2[1]=0.4;
  do
  {
    Shallow_2d_mesh.update( 0.49, std::abs( Shallow_2d_mesh.get_time() - t_end ) );
    Shallow_2d_mesh.dump_gnu( filename_stub + Utility::stringify( file_counter ) + "_gnu.dat" );
    file_counter += 1;
  } while ( Shallow_2d_mesh.get_time() < t_end );

  double h_clawpack( 1.13466 );
  double theta = M_PI / 4;
  DenseVector<double> x( 2, 0.0 ); x[ 0 ] = 0.5 * cos( theta ); x[ 1 ] = 0.5 * sin( theta );
  double h_diag = Shallow_2d_mesh.get_point_values( x )[ h ];
  if ( abs( h_diag - h_clawpack ) > 1.e-3 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " deviation from the Clawpack data = " << abs( h_diag - h_clawpack ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
