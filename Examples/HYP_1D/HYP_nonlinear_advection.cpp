/// \file HYP_nonlinear_advection.cpp
/// \ingroup Examples
/// \ingroup HYP_1D
/// Solving the 1D `nonlinear advection equation'
/// \f[ Q_t + \left ( \frac{Q^2}{2} \right )_x = 0 \quad \mbox{where} \quad Q=Q(x,t) \f] using a TVD Lax-Friedrichs scheme
/// for \f$ x\in[0,1]\f$. The initial condition is taken to be
/// \f[ Q(x,0) =  \sin( 2\pi x ) \f] The test is simply conservation of \f$ Q \f$ in this case.

#include <OneD_HYP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {

    /// Define the system
    class NlinAdv : public OneD_Hyperbolic_System
    {

    public:

      /// One dimemsional inviscid Burgers problem
      NlinAdv() : OneD_Hyperbolic_System( 1 )
      {}

      /// Define the vector flux
      void flux_fn( const double& x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ 0 ] = q[ 0 ] * q[ 0 ] / 2;
      }

      /// Bound the shock speed
      double max_charac_speed( const DenseVector<double> &q ) const
      {
        // maximum shock speed
        return std::abs( q[ 0 ] );
      }

      ///// edge conditions
      std::vector<bool> edge_values( const int face_index, const double& x, DenseVector<double>& q, const double &t  ) const
      {
        std::vector<bool> inflow( q.size(), false );
        // x doesn't matter since the conditions are fixed
        if ( face_index == -1 )
        {
          q[ 0 ] = 0.0;
          inflow[ 0 ] = true;
        }
        if ( face_index == 1 )
        {
          q[ 0 ] = 0.0;
          inflow[ 0 ] = true;
        }
        return inflow;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, DenseVector<double> &q )
    {
      q[ 0 ] = std::sin( 2 * M_PI * x );
    }
  } //end Example namespace

} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== Hyperbolic: 1D nonlinear advection equation  ====\n";
  cout << "\n";

  // define the domain/mesh
  const double left =  0.0;
  const double right = 1.0;
  const unsigned N = 141;
  DenseVector<double> faces_x = Utility::uniform_node_vector( left, right, N );

  double t = 0.0;

  Example::NlinAdv conservative_problem;
  OneD_TVDLF_Mesh NlinAdv_mesh( faces_x, &conservative_problem, Example::Q_init );
  NlinAdv_mesh.set_limiter( 0 );

  double I1 = NlinAdv_mesh.integrate()[0];
  int loop_counter( 5 );
  TrackerFile my_file( "./DATA/HYP_NlinAdv.dat" );
  // first column of the output will always be the time
  my_file.push_ptr( &t, "time" );
  OneD_Node_Mesh<double> soln = NlinAdv_mesh.get_soln();
  my_file.push_ptr( &soln, "mesh" );

  double asym( 0.0 );
  do
  {
    double dt = NlinAdv_mesh.update( 0.475 );
    t += dt;
    soln = NlinAdv_mesh.get_soln();

    if ( loop_counter % 10 == 0 )
    {
      soln = NlinAdv_mesh.get_soln();
      my_file.update();
      my_file.newline();
    }
    ++loop_counter;
    asym = std::max( asym, std::abs( soln.get_interpolated_vars( 0.75 )[0] + soln.get_interpolated_vars( 0.25 )[0] ) );
  }
  while ( ( t < 0.4 ) && ( loop_counter < 1000 ) );

  double I2 = NlinAdv_mesh.integrate()[0];
  soln = NlinAdv_mesh.get_soln();
  my_file.update();
  my_file.newline();
  // problem should be antisymmetric about x = 1/2
  if ( ( asym > 1.e-10 ) || ( std::abs( I1 - I2 ) > 1.e-8 ) || ( loop_counter >= 1000 ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "asymmetry = " << asym << "\n";
    cout << "integral difference = " << I1 - I2 << "\n";
    cout << "loop counter = " << loop_counter << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

} // end of main()
