/// \file HYP_acoustic_reflection.cpp
/// \ingroup Examples
/// \ingroup HYP_1D
/// Solve the constant coefficient linear sound wave problem
/// \f[ p_t + K u_x = 0 \f]
/// \f[ \rho u_t + p_x = 0 \f]
/// for a right-propagating square pulse in a medium
/// with constant bulk modulus \f$K=1\f$ and density \f$\rho = 1\f$
/// and reflecting boundary conditions at both sides.

#include <string>
#include <cassert>

#include <Utility.h>
#include <OneD_HYP_bundle.h>

enum { p, u };

namespace CppNoddy
{
  namespace Example
  {
    /// bulk modulus
    double K( 2.0 );
    /// higher density
    double rho( 2.0 );

    /// Define the system
    class Acoustic_1d_ref : public OneD_Hyperbolic_System
    {

    public:

      /// One dimemsional constant coefft acoustic problem
      Acoustic_1d_ref() : OneD_Hyperbolic_System( 2 )
      {}

      /// Define the vector flux
      void flux_fn( const double&x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ p ] = K * q[ u ];
        f[ u ] = q[ p ] / rho;
      }

      /// Bound the shock speed
      double max_charac_speed( const DenseVector<double> &q ) const
      {
        // sound speeds
        double c = sqrt( K / rho );
        // maximum shock speed
        return c;
      }

      /// edge conditions
      std::vector<bool> edge_values( const int face_index, const double& x, DenseVector<double>& q ) const
      {
        // reflection condition
        std::vector<bool> inflow( q.size(), false );
        q[ u ] = 0.0;
        inflow[ u ] = true;
        return inflow;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      if ( ( x < 0.5 ) && ( x > -0.5 ) )
      {
        q[ p ] = std::exp( -20 * x * x );
        q[ u ] = q[ p ] / sqrt( K * rho );
      }
      else
      {
        q[ p ] = q[ u ] = 0.0;
      }
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
  const double left =  -3.0;
  const double right = 3.0;
  const unsigned N = 800;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  // time
  double t = 0;
  double t_end = 6;

  // hyperbolic problem
  Example::Acoustic_1d_ref conservative_problem;
  OneD_TVDLF_Mesh Acoustic_mesh( faces_x, &conservative_problem, Example::Q_init );
  Acoustic_mesh.set_limiter( 1 );

  // output
  double I1 = Acoustic_mesh.integrate()[0];
  int loop_counter( 0 );
  int file_counter( 1 );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_acoustic_ref";
  TrackerFile my_file( 5 );
  OneD_Node_Mesh<double> soln = Acoustic_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    if ( loop_counter % 10 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Acoustic_mesh.get_soln();
      my_file.update();
      file_counter += 1;
    }
    // take a step
    t += Acoustic_mesh.update( 0.49, t_end - t );
    ++loop_counter;
  }
  while ( ( t < t_end ) && ( loop_counter < 1100 ) );

  double I2 = Acoustic_mesh.integrate()[0];
  // compute the pressure difference after the pressure pulse
  // has returned to its original position
  soln = Acoustic_mesh.get_soln();
  my_file.update();
  double diff = 0.0;
  for ( std::size_t i = 0; i < soln.get_nnodes(); ++i )
  {
    // get the initial condition
    DenseVector<double> q( 2, 0.0 );
    Example::Q_init( soln.coord( i ), q );
    // difference between initial and final state
    diff = std::max( diff, q[0] - soln( i, 0 ) );
  }
  // check the maximum difference in the reflected pulse & the
  // integral of the pressure over the mesh
  if ( ( diff > 1.e-2 ) || ( std::abs( I1 - I2 ) > 1.e-8 ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "difference in the reflected wave = " << diff << "\n";
    cout << "difference in the integral = " << std::abs( I1 - I2 ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

} // end of main()
