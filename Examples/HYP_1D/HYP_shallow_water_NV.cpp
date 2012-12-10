/// \file HYP_shallow_water_NV.cpp
/// \ingroup Examples
/// \ingroup HYP_1D
/// Solve the shallow water equations in one dimension for an
/// initial (small) hump of fluid
/// \f[ h_t +  (uh)_x = 0 \f]
/// \f[ (uh)_t + (hu^2 + gh^2 /2 )_x = 0 \f]
/// where \f$g=1\f$ and the boundaries are reflecting.

#include <string>

#include <Utility.h>
#include <OneD_HYP_bundle.h>

enum { h, hu };

namespace CppNoddy
{
  namespace Example
  {
    /// gravitational acceleration
    double g( 1.0 );
    /// initial hump amplitude
    double A( 0.05 );

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
      std::vector<bool> edge_values( const int face_index, const double& x, DenseVector<double>& q, const double &t  ) const
      {
        std::vector<bool> inflow( q.size(), false );
        q[ hu ] = 0.0;
        inflow[ hu ] = true;
        return inflow;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      q[ h ] = 1 + A * std::exp( -10 * x * x );
      q[ hu ] = 0;
    }
  } //end Example namespace

} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== Hyperbolic: 1D acoustic wave reflection problem =\n";
  cout << "\n";

  // define the domain/mesh
  const double left =  -2.0;
  const double right = 2.0;
  const unsigned N = 512;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  // time
  double t = 0;
  double t_end = 10;

  // hyperbolic problem
  Example::Shallow_1d_ref conservative_problem;
  OneD_TVDLF_Mesh Shallow_mesh( faces_x, &conservative_problem, Example::Q_init );
  Shallow_mesh.set_limiter( 0 );

  // output
  int loop_counter( 4 );
  int file_counter( 1 );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_shallow";
  TrackerFile my_file( 5 );
  OneD_Node_Mesh<double> soln = Shallow_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    t += Shallow_mesh.update( 0.49, t_end - t );

    if ( loop_counter % 5 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Shallow_mesh.get_soln();
      my_file.update();
      file_counter += 1;
    }
    ++loop_counter;
  }
  while ( ( t < t_end ) && ( loop_counter < 3000 ) );

} // end of main()
