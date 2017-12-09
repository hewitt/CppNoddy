/// \file BVP_Berman.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Solving the Berman suction-channel solution in the form
/// \f[ f^(iv)(y) = Re ( f(y)f'''(y) - f'(y)f''(y) ) \f]
/// subject to \f$ f(\pm 1) = \pm 1 \f$ and \f$ f'(\pm 1) = 0 \f$
/// by applying the ODE_BVP class.


#include <BVP_bundle.h>

// enumerate the variables in the ODE
enum {f, fd, fdd, fddd };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Berman equation by inheriting the Equation base class
    class Berman_equation : public Equation_1matrix<double>
    {
    public:

      /// The Berman equation is a 4th order real ODE
      Berman_equation() : Equation_1matrix<double>( 4 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<double> &z, DenseVector<double> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = z[ fdd ];
        g[ fdd ] = z[ fddd ];
        g[ fddd ] = Re * ( z[ f ] * z[ fddd ] - z[ fd ] * z[ fdd ] );
      }

      void matrix0( const DenseVector<double>&x, DenseMatrix<double> &m ) const
      {
        Utility::fill_identity(m);
      }

      // The Reynolds number
      double Re;
    };

    class Berman_left_BC : public Residual<double>
    {
    public:
      // 2 boundary conditions and 4 unknowns
      Berman_left_BC() : Residual<double> ( 2, 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ f ] + 1.0;
        B[ 1 ] = z[ fd ];
      }
    };

    class Berman_right_BC : public Residual<double>
    {
    public:
      // 2 boundary conditions and 4 unknowns
      Berman_right_BC() : Residual<double> ( 2, 4 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &B ) const
      {
        B[ 0 ] = z[ f ] - 1.0;
        B[ 1 ] = z[ fd ];
      }
    };

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: finite-difference solution of Berman eqn ===\n";
  cout << "\n";

  Example::Berman_equation problem;
  Example::Berman_left_BC BC_left;
  Example::Berman_right_BC BC_right;

  // set the Reynolds number
  problem.Re = 1.0;
  // Reynolds number step
  double delta_Re = 0.1;
  // domain is -1 to 1
  double left = -1.0;
  double right = 1.0;
  // number of nodal points
  int N = 1401;
  // mesh
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
  // Example tolerance
  const double tol = 1.e-4;

  // pass it to the ode
  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );

  // our initial guess
  for ( int i = 0; i < N; ++i )
  {
    double y = ode.solution().coord( i );
    // set f(y)
    ode.solution()( i, f ) = 1.5 * ( y - y * y * y / 3 );
    // set f'(y)
    ode.solution()( i, fd ) = 1.5 * ( 1 - y * y );
    // set f''(y)
    ode.solution()( i, fdd ) = -3 * y;
    // set f''(y)
    ode.solution()( i, fddd ) = -3;
  }

  // zero-order continuation in Re
  do
  {
    // solve the problem using 2nd order finite-difference
    try
    {
      try
      {
        ode.solve2();
        // cout << problem.Re << " " << ode.solution().get_interpolated_vars( 0.25 )[0] << "\n";
      }
      catch ( std::runtime_error )
      {
        cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
        return 1;
      }
    }
    catch ( ExceptionBifurcation )
    {
      cout << " Bifurcation detected between Re = " << problem.Re - delta_Re
           << " and Re = " << problem.Re << "\n";
      cout << " Continuing further.\n";
    }
    problem.Re += delta_Re;

  }
  while ( problem.Re < 10.0 + tol );


  // compare against the Shoot_Berman data at Re = 10
  // data to be compared is f''(y) at y = -1
  if ( std::abs( 5.99898 - ode.solution()( 0, fdd ) ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << std::abs( 5.99898 - ode.solution()( 0, fdd ) ) << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
}
