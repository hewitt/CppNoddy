/// \file BVP_JH_arc_NV.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Arc length conitinuation of the equation for self-similar flow
/// between two planar angled walls; the Jeffery-Hamel flow. The code
/// continues the primary branch from Re = 40 around the saddle-node
/// and back to Re = 40.

#include <BVP_bundle.h>

// enumerate the variables
enum { G, Gd, Gdd, Gddd };

namespace CppNoddy
{
  namespace Example
  {

    double alpha;
    double Re;

    class JH_base_equation : public Equation_1matrix<double>
    {
    public:

      JH_base_equation() : Equation_1matrix<double> ( 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[ G ] = z[ Gd ];
        f[ Gd ] = z[ Gdd ];
        f[ Gdd ] = z[ Gddd ];
        f[ Gddd ]  = - 2 * alpha * Re * z[ Gd ] * z[ Gdd ] - 4 * alpha * alpha * z[ Gdd ];
      }
      
      void matrix0( const DenseVector<double>&x, DenseMatrix<double> &m ) const
      {
        Utility::fill_identity(m);
      }
      
    };

    class JH_left_BC : public Residual<double>
    {
    public:
      JH_left_BC() : Residual<double> ( 2, 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ G ] + 0.5;
        B[ 1 ] = z[ Gd ];
      }
    };

    class JH_right_BC : public Residual<double>
    {
    public:
      JH_right_BC() : Residual<double> ( 2, 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ G ] - 0.5;
        B[ 1 ] = z[ Gd ];
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{

  // set up the problem
  Example::JH_base_equation problem;
  Example::alpha = 0.1;
  Example::Re = 5.0;
  // set up BCs
  Example::JH_left_BC BC_left;
  Example::JH_right_BC BC_right;
  // number of points
  unsigned N = 201;
  // set up the domain from -1 to 1
  double left = -1.0;
  double right = 1.0;
  // mesh
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );

  // pass it to the ode
  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );

  // output file
  TrackerFile my_file( "./DATA/BVP_JH.dat" );
  my_file.header();
  double U_centre;
  my_file.push_ptr( &Example::Re, "Re" );
  my_file.push_ptr( &U_centre, "Centerline velocity" );

  // solve the base flow
  ode.solve2();

  // initialise the arc-length routine
  double ds( 0.01 );
  ode.set_monitor_det( false );
  ode.init_arc( &Example::Re, ds, 0.05 );
  ode.rescale_theta() = true;
  ode.desired_arc_proportion() = 0.5;

  do
  {
    ds = ode.arclength_solve( ds );
    U_centre = ode.solution().get_interpolated_vars( 0.0 )[ Gd ];
    my_file.update();
  }
  while ( Example::Re > 5 );

}
