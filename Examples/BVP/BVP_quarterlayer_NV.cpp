/// \file BVP_Harmonic.cpp
/// \ingroup Examples
/// \ingroup BVP
/// Solving the (somewhat ad-hoc) nonlinear version of the Stewatson quarter layer
/// as proposed by Barcilon (1970). Somewhat ad-hoc because it assumes the linear
/// relation for Ekman pumping is still valid despite the finite Rossby number.
/// \f[ V''(\xi) = (V(\xi)-V_\infty ) \left ( 1 + \lambda V'(\xi) ) \f]
/// subject to \f$ V(0) = 0 \f$ and \f$ V(\infty) \to V_\infty \f$.

#include <cassert>

#include <BVP_bundle.h>

// enumerate the variables in the ODE
enum {V, Vd };

namespace CppNoddy
{
  namespace Example
  {
    
    double Vinf=2;
    double lambda=0;//reset below using Ekman number and Rossby number
    
    /// Define the harmonic equation by inheriting the Equation base class
    template <typename _Type, typename _Xtype>
    class Barcilon_equation : public Equation_1matrix<_Type, _Xtype>
    {
    public:

      /// The equation is a 2nd order ODE
      Barcilon_equation() : Equation_1matrix<_Type, _Xtype>( 2 ) {}

      /// The Berman equation
      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &g ) const
      {
        g[ V ] = z[ Vd ];
        g[ Vd ] = ( z[ V ] - Vinf )*( 1 + lambda*z[ Vd ] );
      }
      
      void matrix0( const DenseVector<_Type>&x, DenseMatrix<_Type> &m ) const
      {
        Utility::fill_identity(m);
      }

    };

    template <typename _Type>
    class Barcilon_right_BC : public Residual<_Type>
    {
    public:
      Barcilon_right_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ V ];
      }
    };

    template <typename _Type>
    class Barcilon_left_BC : public Residual<_Type>
    {
    public:
      Barcilon_left_BC<_Type>() : Residual<_Type> ( 1, 2 ) {}

      void residual_fn( const DenseVector<_Type> &z, DenseVector<_Type> &B ) const
      {
        B[ 0 ] = z[ V ] - Vinf;
      }
    };


  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP: finite-difference solution of Barcilon's nonlinear quarter layer eqn  =\n";
  cout << "\n";

  Example::Barcilon_equation<double, double> problem;
  Example::Barcilon_left_BC <double> BC_left;
  Example::Barcilon_right_BC<double> BC_right;

  double left =  -30.0;
  double right =  0.0;
  // number of nodal points
  unsigned N = 201;
  // a real mesh -- along real axis
  DenseVector<double> xi_nodes( Utility::uniform_node_vector( left, right, N ) );

  // a real ODE BVP
  ODE_BVP<double> ode( &problem, xi_nodes, &BC_left, &BC_right );

  // our initial guess
  for ( unsigned i = 0; i < N; ++i )
  {
    double xi = ode.solution().coord( i );
    // set V(xi)
    ode.solution()( i, V ) = Example::Vinf*(1+exp(xi));
    // set V'(xi)
    ode.solution()( i, Vd ) = Example::Vinf*exp(xi);
  }


  double Ek( 1.e-4 );
  double Ro( 0.0 );
  TrackerFile metric("./DATA/Barcilon.dat");
  metric.push_ptr( &Ro, "Ro" );
  metric.push_ptr( &Example::lambda, "lambda" );
  metric.header();
  do
  { 
    Example::lambda = Ro*pow( sqrt(Ek)*Example::Vinf/2., -0.5 );     
    ode.solve2();  
    metric.update();
    ode.solution().dump_gnu("./DATA/Barcilon_" + Utility::stringify( Ro, 4 ) + ".dat");
    cout << "lambda = " << Example::lambda << " Ek = " << Ek << " Ro = " << Ro << "\n";
    Ro += 0.005;  
  } while ( Ro < 1.0 );
}
