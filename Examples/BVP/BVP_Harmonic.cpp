/// \file BVP_Harmonic.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Solving the Harmonic equation
/// \f[ f''(z) + f(z) = 0 \f]
/// subject to \f$ f(0) = 0 \f$ and \f$ f(1) = 1 \f$ (OR \f$ f(i) = 1 \f$)
/// by applying the ODE_BVP class. We solve for two different cases, the
/// first is a real BVP along the real axis, the second is a complex BVP
/// along the imaginary axis.

#include <cassert>

#include <BVP_bundle.h>

// enumerate the variables in the ODE
enum {f, fd };

namespace CppNoddy
{
  namespace Example
  {
    /// Define the harmonic equation by inheriting the Equation base class
    template <typename _Type, typename _Xtype>
    class Harmonic_equation : public Equation_1matrix<_Type, _Xtype>
    {
    public:

      /// The harmonic equation is a 2nd order ODE
      Harmonic_equation() : Equation_1matrix<_Type, _Xtype>( 2 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = - z[ f ];
      }
      
      void matrix0( const DenseVector<_Type>&x, DenseMatrix<_Type> &m ) const
      {
        Utility::fill_identity(m);
      }

    };

    template <typename _Type>
    class Harmonic_left_BC : public Residual<_Type>
    {
    public:
      Harmonic_left_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ f ];
      }
    };

    template <typename _Type>
    class Harmonic_right_BC : public Residual<_Type>
    {
    public:
      Harmonic_right_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ f ] - 1.;
      }
    };


  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: finite-difference solution of Harmonic eqn =\n";
  cout << "\n";

  Example::Harmonic_equation<double, double> real_problem;
  Example::Harmonic_left_BC <double> real_BC_left;
  Example::Harmonic_right_BC<double> real_BC_right;
  Example::Harmonic_equation<D_complex, D_complex> complex_problem;
  Example::Harmonic_left_BC <D_complex> complex_BC_left;
  Example::Harmonic_right_BC<D_complex> complex_BC_right;

  double left =  0.0;
  double right = 1.0;
  // number of nodal points
  unsigned N = 21;
  // a real mesh -- along real axis
  DenseVector<double> real_nodes( Utility::uniform_node_vector( left, right, N ) );
  // a complex mesh -- along imaginary axis
  D_complex eye( 0.0, 1.0 );
  DenseVector<D_complex> complex_nodes( real_nodes );
  complex_nodes *= eye;
  // Example tolerance
  const double tol = 1.e-4;

  // a real ODE BVP
  ODE_BVP<double> real_ode( &real_problem, real_nodes, &real_BC_left, &real_BC_right );
  // a complex ODE BVP solved on a line in the complex plane
  ODE_BVP<D_complex, D_complex> complex_ode( &complex_problem, complex_nodes, &complex_BC_left, &complex_BC_right );
  complex_ode.set_monitor_det( false );

  // our initial guess
  for ( unsigned i = 0; i < N; ++i )
  {
    double y = real_ode.solution().coord( i );
    // set f(y)
    real_ode.solution()( i, f ) = y;
    complex_ode.solution()( i, f ) = y;
    // set f'(y)
    real_ode.solution()( i, fd ) = 1;
    complex_ode.solution()( i, fd ) = 1;
  }

  // solve the problem using 2nd order finite-difference
  try
  {
    real_ode.solve2();
    complex_ode.solve2();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  // find the maximum deviation from the analytical solutions of sin(x)/sin(1) and sinh(z)/sinh(1)
  double real_diff = 0;
  double complex_diff = 0;
  for ( unsigned i = 0; i < N; ++i )
  {
    real_diff = std::max( std::abs( real_ode.solution()( i, f ) - sin( real_ode.solution().coord( i ) ) / sin( 1 ) ), real_diff );
    complex_diff = std::max( std::abs( complex_ode.solution()( i, f ) - sin( complex_ode.solution().coord( i ) ) / sin( eye ) ), complex_diff );
  }

  // validation test
  if ( real_diff > tol || complex_diff > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "Real problem error " << real_diff << "\n";
    cout << "Complex problem error " << complex_diff << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
