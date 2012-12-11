/// \file BVP_Troesch.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Solving the Harmonic equation
/// \f[ f''(y) = c \sinh (cy)  \f]
/// subject to \f$ f(0) = 0 \f$ and \f$ f(1) = 1 \f$.

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
    class Troesch_equation : public Equation_1matrix<_Type, _Xtype>
    {
    public:
      double c;

      /// The harmonic equation is a 2nd order ODE
      Troesch_equation() : Equation_1matrix<_Type, _Xtype>( 2 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &g ) const
      {
        g[ f ] = z[ fd ];
        g[ fd ] = c*std::sinh(c*z[f]);
      }
      
      void matrix0( const DenseVector<_Type>&x, DenseMatrix<_Type> &m ) const
      {
        Utility::fill_identity(m);
      }

    };

    template <typename _Type>
    class Troesch_left_BC : public Residual<_Type>
    {
    public:
      Troesch_left_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ f ];
      }
    };

    template <typename _Type>
    class Troesch_right_BC : public Residual<_Type>
    {
    public:
      Troesch_right_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

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
  cout << "=== BVP: finite-difference solution of Troesch eqn =\n";
  cout << "\n";

  Example::Troesch_equation<double, double> problem;
  Example::Troesch_left_BC <double> BC_left;
  Example::Troesch_right_BC<double> BC_right;

  double left =  0.0;
  double right = 1.0;
  // number of nodal points
  unsigned N = 151;
  // a real mesh -- along real axis
  DenseVector<double> nodes( Utility::uniform_node_vector( left, right, N ) );
  // Example tolerance
  const double tol = 1.e-4;

  // a real ODE BVP
  problem.c = 5.0;
  ODE_BVP<double> ode( &problem, nodes, &BC_left, &BC_right );

  // our initial guess
  for ( unsigned i = 0; i < N; ++i )
  {
    double y = ode.solution().coord( i );
    // set f(y)
    ode.solution()( i, f ) = y;
    // set f'(y)
    ode.solution()( i, fd ) = 1;
  }

  // solve the problem using 2nd order finite-difference
  try
  {
    ode.solve2();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  //ode.solution().dump_gnu("./DATA/test.dat");

  // find deviation from y=0.5 solution obtianed from: http://www2.imperial.ac.uk/~jcash/BVP_software/readme.php
  const double diff = std::abs( ode.solution().get_interpolated_vars(0.5)[f] - 0.55437396E-01 );

  // validation test
  if ( diff > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "Real problem error " << diff << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
