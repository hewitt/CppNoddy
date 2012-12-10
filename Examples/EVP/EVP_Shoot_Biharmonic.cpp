/// \file EVP_Shoot_Biharmonic.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solving a one-dimensional "Bi-harmonic" eigenvalue problem (EVP)
///  \f[ \left ( \frac{\mbox{d}^2}{\mbox{d}x^2} - \lambda \right )^2 f(x) = 0\,, \quad \mbox{where} \quad f(0)=f(1)=f'(0)=f'(1)=0\,,\f]
/// via Runge-Kutta and (vector) Newton iteration. The problem is
/// nonlinear since \f$\lambda\f$ is unknown.
///
/// Locally solving for eigenvalues using IVP methods is generally speaking
/// not a great way of doing things; this is just a Example rather than
/// any recommendation! See the other EVP examples for better methods.

#include <cassert>

#include <IVP_bundle.h>
#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the biharmonic eigenvalue equation by
    /// inheriting Equation base class.
    class Biharmonic_equation : public Equation<D_complex>
    {
    public:
      /// The biharmonic is a 4th order complex problem
      Biharmonic_equation() : Equation<D_complex>( 4 ) {}

      /// We implement the equation as 4 first-order ODEs.
      void residual_fn( const DenseVector<D_complex>& z, DenseVector<D_complex>& f ) const
      {
        f[ 0 ] = z[ 1 ];
        f[ 1 ] = z[ 2 ];
        f[ 2 ] = z[ 3 ];
        f[ 3 ] = - z[ 2 ] * 2.0 * lambda - z[ 0 ] * pow( lambda, 2 );
      }

      /// The eigenvalue
      D_complex lambda;

    };

    /// Define a residual function using the boundary
    /// conditions for the biharmonic_equation.
    class Biharmonic_residual : public Residual<D_complex>
    {
    public:
      ODE_IVP<D_complex> *ode;

      Biharmonic_residual() : Residual<D_complex> ( 2 )
      {
        // instantiate the 1D bi-harmonic eqn
        eqn = new Biharmonic_equation;
        // set up the ODE domain, step etc
        ode = new ODE_IVP<D_complex>( eqn, 0.0, 1.0, 1000 );
      }

      ~Biharmonic_residual()
      {
        delete ode;
        delete eqn;
      }

      /// implement a residual function.
      void residual_fn( const DenseVector<D_complex>& unknown, DenseVector<D_complex>& BC ) const
      {
        DenseVector<D_complex> u( 4, 0.0 );        // 4th order problem
        u[ 0 ] = 0.0;                 // f = 0
        u[ 1 ] = 0.0;                 // f' = 0
        u[ 2 ] = 1.0;                 // f'' = 1
        u[ 3 ] = unknown[ 0 ];        // f''' = unknown[0]
        eqn -> lambda = unknown[ 1 ]; // parameter = unknown[1]
        // shoot using R-K fixed step size algorithm
        u = ode -> shoot4( u );
        // return the BC residuals
        BC[ 0 ] = u[ 0 ];
        BC[ 1 ] = u[ 1 ]; // trying to make f = f' = 0 at other boundary
      }
    private:
      Biharmonic_equation *eqn;
    };
  } //end Example namespace
} //end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== EVP: Biharmonic equation via shooting  ==========\n";
  cout << "\n";

  Example::Biharmonic_residual problem;
  // instantiate a complex Vector Newton class for the Biharmonic
  Newton<D_complex> Biharm( &problem );

  const D_complex A( -3.0, 7.0 );        // A guess of f'''(0)
  const D_complex B( 2.21, 1.25 );       // A guess of the eigenvalue
  // rig up an initial guess
  DenseVector<D_complex> unknowns( 2, 0.0 );          // Two complex quantities to to find
  unknowns[ 0 ] = A;
  unknowns[ 1 ] = B * B;

  try
  {
    Biharm.iterate( unknowns );
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  D_complex eigenvalue( std::sqrt( unknowns[ 1 ] ) );

  // output for the data
  TrackerFile my_file( "./DATA/Biharmonic.dat" );
  my_file.push_ptr( &eigenvalue, "eigenvalue" );
  my_file.push_ptr( &problem.ode -> get_mesh(), "biharmonic" );
  my_file.header();
  my_file.update();

  const double tol = 1.e-8; // A pass/fail tolerance
  if ( abs( sin( eigenvalue ) + eigenvalue ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout << "   Diff = " << abs( sin( eigenvalue ) + eigenvalue ) << "\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }

}
