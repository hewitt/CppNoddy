/// \file BVPNonIdentityMatrix0.cpp
/// \ingroup Tests
/// \ingroup BVP
/// Solving the equation
/// \f[ f(y) f''(y) + f'(y)^2 = 1+\gamma y  \f]
/// subject to \f$ f(0) = 1 \f$ and \f$ f(1) = 2 \f$.
/// We don't divide by the \f$ f(y) \f$ that multiplies the
/// highest derivative, rather we define a non-identity matrix0.
/// The solution is compared at all nodes to the exact solution
/// \f[ f(y)=\left ( y^2 + \frac{\gamma}{3}y^3 + (2-\frac{\gamma}{3}) y + 1 \right )^{1/2}\,.\f]

#include <BVP_bundle.h>

// enumerate the variables in the ODE
enum {f, fd };

namespace CppNoddy
{
  namespace Example
  {

    /// Define the harmonic equation by inheriting the Equation base class
    template <typename _Type, typename _Xtype>
    class Nonidentity_equation : public Equation_1matrix<_Type, _Xtype>
    {
    public:

      double gamma;

      /// The harmonic equation is a 2nd order ODE
      Nonidentity_equation() : Equation_1matrix<_Type, _Xtype>( 2 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &g ) const
      {
        const double y( this -> coord(0) );
        g[ f ] = z[ fd ];
        g[ fd ] = 1 + gamma*y - z[ fd ]*z[ fd ];
      }

      void matrix0( const DenseVector<_Type>&z, DenseMatrix<_Type> &m ) const
      {
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = z[f];
      }

    };

    template <typename _Type>
    class Nonidentity_left_BC : public Residual<_Type>
    {
    public:
      Nonidentity_left_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ f ] - 1;
      }
    };

    template <typename _Type>
    class Nonidentity_right_BC : public Residual<_Type>
    {
    public:
      Nonidentity_right_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ f ] - 2;
      }
    };


  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: finite-difference, non-identity matrix0 ====\n";
  cout << "\n";

  Example::Nonidentity_equation<double, double> problem;
  Example::Nonidentity_left_BC <double> BC_left;
  Example::Nonidentity_right_BC<double> BC_right;

  double left =  0.0;
  double right = 1.0;
  // number of nodal points
  unsigned N = 151;
  // a real mesh -- along real axis
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
  // Example tolerance
  const double tol = 1.e-4;

  // a real ODE BVP
  problem.gamma = 5.0;
  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );

  // our initial guess
  for ( unsigned i = 0; i < N; ++i )
  {
    double y = ode.solution().coord( i );
    // set f(y)
    ode.solution()( i, f ) = 1+y;
    // set f'(y)
    ode.solution()( i, fd ) = 1;
  }

  // solve the problem using 2nd order finite-difference
  try
  {
    ode.solve2();
  }
  catch (const std::runtime_error& error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }


  // find deviation from exact solution
  double diff = 0;
  for ( unsigned i = 0; i < N; ++i )
  {
    const double y( ode.solution().coord(i) );
    const double exact( sqrt(y*y + problem.gamma*y*y*y/3 + (2-problem.gamma/3)*y + 1) );
    diff = std::max( std::abs( ode.solution().get_interpolated_vars(y)[f] - exact ), diff );
  }

  // validation test
  if ( diff > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "Real problem error " << diff << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
}
