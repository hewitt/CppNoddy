/// \file HYP_shallow_source_NV.cpp
/// \ingroup Examples
/// \ingroup HYP_1D
/// Solve the shallow water equations in one dimension
/// \f[ h_t +  (uh)_x = 0 \f]
/// \f[ (uh)_t + (hu^2 + gh^2 /2 )_x = -ghz'(x) \f]
/// where the momentum is fixed upstream, and the
/// flow depth is fixed downstream and \f$ z(x) \f$ is
/// the topography shape. A shock in fluid height should
/// occur on the downstream side of the topography.


#include <string>
#include <cassert>

#include <Utility.h>
#include <OneD_HYP_bundle.h>

enum { h, hu };

namespace CppNoddy
{
  namespace Example
  {
    /// gravitational acceleration
    double g( 9.81 );
    /// hump amplitude
    double A( 0.2 );
    /// Topography shape
    double z( const double &x )
    {
      return A * std::exp( -.5 * (x-10) * (x-10) );
    }
    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      q[ h ] = 0.33;
      q[ hu ] = 0.18;
    }

    /// Define the system
    class Shallow_1d_ref : public OneD_Hyperbolic_System
    {

    public:

      /// One dimemsional constant coefft acoustic problem
      Shallow_1d_ref() : OneD_Hyperbolic_System( 2 )
      {}

      /// Define the vector flux
      void flux_fn( const double &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ h ] = q[ hu ];
        f[ hu ] = q[ hu ] * q[ hu ] / q[ h ] + g * q[ h ] * q[ h ] / 2;
      }

      /// Bound the shock speed
      double max_charac_speed( const DenseVector<double> &q ) const
      {
        // sound speeds
        double c = sqrt( g * q[ h ] );
        // flow speed
        double u = q[ hu ] / q[ h ];
        // maximum shock speed
        return std::max( std::abs( u + c ), std::abs( u - c ) );
      }

      /// edge conditions
      std::vector<bool> edge_values( const int face_index, const double& x, DenseVector<double>& q ) const
      {
        std::vector<bool> inflow( q.size(), false );
        if ( face_index < 0 )
        {
          q[ hu ] = 0.18;
          inflow[ hu ] = true;
        }
        else
        {
          q[ h ] = 0.33;
          inflow[ h ] = true;
        }
        return inflow;
      }

      void source_fn( const double& x, const DenseVector<double>& q, const DenseVector<double>& slope, DenseVector<double>& r ) const
      {
        const double delta( 1.e-6 );
        r[ h ] = 0.0;
        r[ hu ] = - g * q[ h ] * ( z( x + delta ) - z( x ) ) / delta;
      }

    };


  } //end Example namespace

} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== Hyperbolic: 1D shallow water solver + sources ===\n";
  cout << "\n";

  // define the domain/mesh
  const double left = 0.0;
  const double right = 25.0;
  const unsigned N = 401;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  // time
  double t = 0;
  double t_end = 5.;

  // hyperbolic problem
  Example::Shallow_1d_ref conservative_problem;
  OneD_TVDLF_Mesh Shallow_mesh( faces_x, &conservative_problem, Example::Q_init );
  Shallow_mesh.set_limiter( 0 );

  // output
  int loop_counter( 0 );
  int file_counter( 0 );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_shallow";
  TrackerFile my_file( 8 );
  OneD_Node_Mesh<double> soln = Shallow_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    if ( loop_counter % 10 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Shallow_mesh.get_soln();
      for ( unsigned i = 0; i < soln.get_nnodes(); ++i )
      {
        soln( i, h ) += Example::z( soln.coord( i ) );
      }
      my_file.update();
      file_counter += 1;
    }
    t += Shallow_mesh.update( 0.49, t_end - t );
    //cout << t << ", h(10) = " << soln.get_interpolated_vars( 10.0 )[ h ] - Example::z( 10.0 ) << "\n";
    ++loop_counter;
  }
  while ( ( t < t_end ) );

} // end of main()
