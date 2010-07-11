/// \file HYP_acoustic_impedance.cpp
/// \ingroup Examples
/// \ingroup HYP_1D
/// Solve the constant coefficient linear sound wave problem
/// \f[ p_t + K u_x = 0 \f]
/// \f[ \rho u_t + p_x = 0 \f]
/// for a right-propagating square pressure pulse in a medium
/// with constant bulk modulus \f$K=4\f$ and density \f$\rho = 1\f$
/// if \f$\vert x \vert < 1\f$ and \f$\rho = 4\f$ elsewhere.
/// The (first) reflected wave's velocity amplitude is
/// \f[ - \left ( \frac{\frac{K_1}{K_2}-\frac{c_1}{c_2}}{\frac{K_1}{K_2}+\frac{c_1}{c_2}} \right ) \f]
/// where \f$ K_{1,2} \f$ and \f$c_{1,2}\f$ are the bulk modulus and wave speed
/// before and after the density interface, with \f$c_i = \sqrt{ K_i \rho_i } \f$.

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
    double K( 4.0 );
    /// higher density
    double rho_big( 4.0 );
    /// lower density
    double rho_small( 1.0 );
    /// Density function for the medium
    double rho( const double &x )
    {
      if ( std::abs( x ) < 1.0 )
      {
        return rho_small;
      }
      else
      {
        return rho_big;
      }
    }

    /// Define the system
    class Acoustic_1d : public OneD_Hyperbolic_System
    {

    public:

      /// One dimemsional constant coefft acoustic problem
      Acoustic_1d() : OneD_Hyperbolic_System( 2 )
      {}

      /// Define the vector flux
      void flux_fn( const double &x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ p ] = K * q[ u ];
        f[ u ] = q[ p ] / rho( x );
      }

      /// Bound the shock speed
      double max_charac_speed( const DenseVector<double> &q ) const
      {
        // sound speeds
        double c = sqrt( K / rho_small );
        // maximum shock speed
        return c;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      if ( ( x < -1.5 ) && ( x > -2.5 ) )
      {
        q[ p ] = 1.0;
        q[ u ] = q[ p ] / sqrt( K * rho( x ) );
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
  cout << "=== Hyperbolic: 1D acoustic wave, impedance problem =\n";
  cout << "\n";

  // define the domain/mesh
  const double left =  -7.0;
  const double right = 7.0;
  const unsigned N = 800;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  // time
  double t = 0.0;
  double t_end = 5.0;

  // hyperbolic problem
  Example::Acoustic_1d conservative_problem;
  OneD_TVDLF_Mesh Acoustic_mesh( faces_x, &conservative_problem, Example::Q_init );
  Acoustic_mesh.set_limiter( 1 );

  // mesh & test info
  double I1 = Acoustic_mesh.integrate()[0];
  int loop_counter( 4 );
  int file_counter( 1 );
  std::string filename_stub;
  filename_stub = "./DATA/HYP_acoustic_imp";
  TrackerFile my_file( 5 );
  OneD_Node_Mesh<double> soln = Acoustic_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  do
  {
    t += Acoustic_mesh.update( 0.499, t_end - t );

    if ( loop_counter % 2 == 0 )
    {
      my_file.set_filename( filename_stub + Utility::stringify( file_counter ) + ".dat" );
      soln = Acoustic_mesh.get_soln();
      my_file.update();
      file_counter += 1;
    }
    ++loop_counter;

  }
  while ( ( t < t_end ) && ( loop_counter < 1500 ) );

  soln = Acoustic_mesh.get_soln();
  double I2 = Acoustic_mesh.integrate()[0];
  // the first reflected and first transmitted pulses can be determined explicitly.
  // they have amplitude ratios of 1/3 and 8/9ths of the input amplitude for the pressure.
  soln = Acoustic_mesh.get_soln();
  double E1 = std::abs( soln.get_interpolated_vars( -5.0 )[0] - 1. / 3. );
  double E2 = std::abs( soln.get_interpolated_vars(  4.0 )[0] - 8. / 9. );
  if ( ( E1 > 1.e-4 ) || ( E2 > 1.e-4 ) || ( std::abs( I1 - I2 ) > 1.e-8 ) || ( loop_counter >= 1500 ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

} // end of main()
