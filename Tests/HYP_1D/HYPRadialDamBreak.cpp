/// \file HYPRadialDamBreak.cpp
/// \ingroup Tests
/// \ingroup HYP_1D
/// Solve the shallow water equations in one dimension for an
/// initial column of fluid
/// \f[ h_t +  (uh)_r = -hu/r \f]
/// \f[ (uh)_t + (hu^2 + gh^2 /2 )_r = -hu^2/r \f]
/// The result is compared to the same problem solved using
/// Clawpack, evaluated at the single point x=0.5.

#include <OneD_HYP_bundle.h>

enum { h, hu };

namespace CppNoddy
{
  namespace Example
  {
    /// gravitational acceleration
    double g( 1.0 );
    /// initial hump amplitude
    double A( 1.0 );

    /// Define the system
    class Shallow_1d_rad : public OneD_Hyperbolic_System
    {

    public:

      /// One dimemsional constant coefft acoustic problem
      Shallow_1d_rad() : OneD_Hyperbolic_System( 2 )
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
        const double c( sqrt( g * q[ h ] ) );
        // flow speed
        const double u( q[ hu ] / q[ h ] );
        // maximum shock speed
        return std::max( std::abs( u + c ), std::abs( u - c ) );
      }

      void source_fn( const double& x, const DenseVector<double>& q, const DenseVector<double>& slope, DenseVector<double>& r ) const
      {
        r[ h ] = -q[ hu ] / x;
        r[ hu ] = -q[ hu ] * q[ hu ] / ( q[ h ] * x );
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      double xi( x - 0.5 );
      q[ h ] = 1 + 0.5 * ( 1 - std::tanh( 100 * xi ) );
      q[ hu ] = 0;
    }
  } //end Example namespace

} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== Hyperbolic: Shallow water radial dam break ======\n";
  cout << "\n";

  // There is a singular source_fn at r=0, so we bodge it here (for now).
  /// \todo Include a mechanism for avoiding source term computations
  /// at the singular point by indicating where edge conditions are specified.
  /// For the time being, we'll bodge it by stopping away from x=0.
  const double left =  1.e-4;
  const double right = 2.5;
  const unsigned N = 2001;
  DenseVector<double> faces_x = Utility::power_node_vector( left, right, N, 1.0 );

  // time
  double t = 0;
  double t_end = 0.5;

  // hyperbolic problem
  Example::Shallow_1d_rad conservative_problem;
  OneD_TVDLF_Mesh Shallow_mesh( faces_x, &conservative_problem, Example::Q_init );
  Shallow_mesh.set_limiter( 0 );

  // output
  int loop_counter( 0 );
  int file_counter( 0 );

  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_shallow_rad";
  TrackerFile my_file( 5 );
  OneD_Node_Mesh<double> soln = Shallow_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    if ( loop_counter % 50 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Shallow_mesh.get_soln();
      my_file.update();
      file_counter += 1;
    }
    t += Shallow_mesh.update( 0.499, t_end - t );
    ++loop_counter;
  }
  while ( ( t < t_end ) && ( loop_counter < 3000 ) );

  // we test the result against the Clawpack computational result
  soln = Shallow_mesh.get_soln();
  double h_clawpack( 1.13466 );
  double h_test( soln.get_interpolated_vars( 0.5 )[ h ] );
  if ( abs( h_test - h_clawpack ) > 1.e-3 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << " deviation from the Clawpack data = " << abs( h_test - h_clawpack ) << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }

} // end of main()
