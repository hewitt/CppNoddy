/// \file HYP_2D_shallow_source.cpp
/// \ingroup Examples
/// \ingroup HYP_2D
/// Two dimensional shallow water equations over topography. In this
/// case, although the computation is 2D the topography is independent
/// of the transverse coordinate and the result is checked against the
/// 1D code.

#include <string>

#include <Utility.h>
#include <TwoD_HYP_bundle.h>

enum { h, hu, hv };

namespace CppNoddy
{
  namespace Example
  {
    /// gravitational acceleration
    double g( 9.81 );
    /// hump amplitude
    double A( 0.2 );
    /// Topography shape
    double z( const double &x, const double &y )
    {
      // to compare with our 1D formulation, we'll make it independent
      // of the transverse y-coordinate
      double rsq( ( x - 10 ) * ( x - 10 ) + 0.0 * y * y );
      return A * std::exp( - 0.5 * rsq );
    }

    /// Set the initial state of the system
    void Q_init( const double &x, const double &y, DenseVector<double> &q )
    {
      q[ h ]  = 0.33;
      q[ hu ] = 0.18;
      q[ hv ] = 0.0;
    }

    /// Define the system
    class Shallow_2d_source : public TwoD_Hyperbolic_System
    {

    public:

      /// One dimemsional constant coefft acoustic problem
      Shallow_2d_source() : TwoD_Hyperbolic_System( 3 )
      {}

      /// Define the vector flux
      void flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ h ] = q[ hu ];
        f[ hu ] = q[ hu ] * q[ hu ] / q[ h ] + g * q[ h ] * q[ h ] / 2;
        f[ hv ] = q[ hu ] * q[ hv ] / q[ h ];
      }

      void flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ h ] = q[ hv ];
        f[ hu ] = q[ hu ] * q[ hv ] / q[ h ];
        f[ hv ] = q[ hv ] * q[ hv ] / q[ h ] + g * q[ h ] * q[ h ] / 2;
      }

      // we'll use analytic Jacobians
      void Jac_flux_fn_x( const DenseVector<double> &x, const DenseVector<double> &q, DenseMatrix<double> &J ) const
      {
        J( 0, hu ) = 1;
        //
        J( 1, h  ) = - q[ hu ] * q[ hu ] / ( q[ h ] * q[ h ] ) + g * q[ h ];
        J( 1, hu ) = 2 * q[ hu ] / q[ h ];
        //
        J( 2, h  ) = - q[ hu ] * q[ hv ] / ( q[ h ] * q[ h ] );
        J( 2, hu ) = q[ hv ] / q[ h ];
        J( 2, hv ) = q[ hu ] / q[ h ];
      }

      void Jac_flux_fn_y( const DenseVector<double> &x, const DenseVector<double> &q, DenseMatrix<double> &J ) const
      {
        J( 0, hv ) = 1;
        //
        J( 1, h  ) = - q[ hu ] * q[ hv ] / ( q[ h ] * q[ h ] );
        J( 1, hu ) = q[ hv ] / q[ h ];
        J( 1, hv ) = q[ hu ] / q[ h ];
        //
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

      /// edge conditions
      std::vector<bool> edge_values( const int& face_index, const DenseVector<double>& x, DenseVector<double>& q ) const
      {
        std::vector<bool> comps( q.size(), false );
        switch ( face_index )
        {
        case 1:
          q[ h ] = 0.33;
          comps[ h ] = true;
          break;
        case 3:
          q[ hu ] = 0.18;
          comps[ hu ] = true;
          break;
        }
        return comps;
      }

      void source_fn( const DenseVector<double>& x, const DenseVector<double>& q, DenseVector<double>& r ) const
      {
        const double delta( 1.e-6 );
        r[ h ] = 0.0;
        r[ hu ] = - g * q[ h ] * ( z( x[0] + delta, x[1] ) - z( x[0], x[1] ) ) / delta;
        r[ hv ] = - g * q[ h ] * ( z( x[0], x[1] + delta ) - z( x[0], x[1] ) ) / delta;
      }

    };
  } //end Example namespace
} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Hyperbolic: 2D shallow water over topography ====\n";
  cout << "\n";

  std::string filename_stub( "./DATA/HYP_2D_shallow_source" );
  // define the domain/mesh
  const double west =  0;
  const double east = 25;
  const double south = -3;
  const double north = 3;
  const unsigned Nx = 401;
  const unsigned Ny = 3;
  DenseVector<double> faces_x = Utility::power_node_vector( west, east, Nx, 1.0 );
  DenseVector<double> faces_y = Utility::power_node_vector( south, north, Ny, 1.0 );

  Example::Shallow_2d_source conservative_problem;
  TwoD_TVDLF_Mesh Shallow_2d_mesh( faces_x, faces_y, &conservative_problem, Example::Q_init );
  Shallow_2d_mesh.set_limiter( 1 );

  int file_counter( 0 );
  int loop_counter( 0 );
  const double t_end = 5;
  do
  {
    if ( loop_counter % 10 == 0 )
    {
      Shallow_2d_mesh.dump_gnu( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      file_counter += 1;
    }
    Shallow_2d_mesh.update( 0.49, std::abs( Shallow_2d_mesh.get_time() - t_end ) );
    loop_counter += 1;

  }
  while ( Shallow_2d_mesh.get_time() < t_end );

  DenseVector<double> x( 2, 0.0 );
  x[ 0 ] = 10;
  x[ 1 ] = 0;
  double h_10( Shallow_2d_mesh.get_point_values( x )[0] );
  // compare this value to the 1D code's result for this resolution and CFL
  double E( abs( h_10 - 0.116773 ) );
  if ( E > 1.e-5 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << E << " \n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
