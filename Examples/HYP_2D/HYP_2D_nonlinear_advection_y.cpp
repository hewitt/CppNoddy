/// \file HYP_2D_nonlinear_advection_y.cpp
/// \ingroup Examples
/// \ingroup HYP_2D
/// Solving the 1D `nonlinear advection equation'
/// \f[ Q_t + \left ( \frac{Q^2}{2} \right )_y = 0 \quad \mbox{where} \quad Q=Q(x,y,t) \f]
/// using a TVD Lax-Friedrichs scheme
/// for \f$ (x,y)\in[-1,1]\times[-1,1]\f$. The initial condition is a sine distribution.

#include <TwoD_HYP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {

    /// Define the system
    class NlinAdv : public TwoD_Hyperbolic_System
    {

    public:

      /// One dimemsional inviscid Burgers problem
      NlinAdv() : TwoD_Hyperbolic_System( 1 )
      {}

      void flux_fn_x( const DenseVector<double>& x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ 0 ] = 0;
      }

      void flux_fn_y( const DenseVector<double>& x, const DenseVector<double> &q, DenseVector<double> &f ) const
      {
        f[ 0 ] = q[ 0 ] * q[ 0 ] / 2;
      }


      /// Bound the wave speed
      void max_charac_speed( const DenseVector<double> &x, const DenseVector<double> &q, DenseVector<double> &c ) const
      {
        // maximum wave speed
        c[ 0 ] = c[ 1 ] = std::abs( q[ 0 ] );
      }

      ///// edge conditions
      std::vector<bool> edge_values( const int& face_index, const DenseVector<double>& x, DenseVector<double>& q ) const
      {
        std::vector<bool> inflow( q.size(), false );
        // x doesn't matter since the conditions are fixed
        if ( face_index == 1 )
        {
          q[ 0 ] = 0.0;
          inflow[ 0 ] = true;
        }
        if ( face_index == 3 )
        {
          q[ 0 ] = 0.0;
          inflow[ 0 ] = true;
        }
        return inflow;
      }

    };

    /// Set the initial state of the system
    void Q_init( const double &x, const double &y, DenseVector<double> &q )
    {
      q[ 0 ] = std::sin( 2 * M_PI * y );
    }
  } //end Example namespace

} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Hyperbolic: 2D nonlinear advection in y =========\n";
  cout << "\n";

  // define the domain/mesh
  const double west =  1.0;
  const double east = 0.0;
  const double south =  0.0;
  const double north = 1.0;
  const unsigned N = 51;
  DenseVector<double> faces_x = Utility::uniform_node_vector( east, west, N );
  DenseVector<double> faces_y = Utility::uniform_node_vector( south, north, N );

  Example::NlinAdv conservative_problem;
  TwoD_TVDLF_Mesh NlinAdv_mesh( faces_x, faces_y, &conservative_problem, Example::Q_init );
  NlinAdv_mesh.set_limiter( 0 );

  double asym( 0.0 );
  unsigned loop_counter( 0 );
  DenseVector<double> x1( 2, 0.0 );
  x1[ 0 ] = 0.5;
  x1[ 1 ] = 0.75;
  DenseVector<double> x2( 2, 0.0 );
  x2[ 0 ] = 0.5;
  x2[ 1 ] = 0.25;
  do
  {
    NlinAdv_mesh.update( 0.49 );
    asym = std::max( asym, std::abs( NlinAdv_mesh.get_point_values( x1 )[0] + NlinAdv_mesh.get_point_values( x2 )[0] ) );
    ++loop_counter;
  }
  while ( ( NlinAdv_mesh.get_time() < 0.4 ) && ( loop_counter < 1000 ) );

  // problem should be antisymmetric about y = 1/2
  if ( ( asym > 1.e-9 ) || ( loop_counter >= 1000 ) )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "asymmetry = " << asym << "\n";
    cout << "loop counter = " << loop_counter << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

} // end of main()
