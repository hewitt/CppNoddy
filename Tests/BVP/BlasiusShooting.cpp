/// \file BlasiusShooting.cpp
/// \ingroup Tests
/// \ingroup BVP
/// Solving the Blasius equation
/// \f[ f'''(y) + f(y) f''(y) = 0\,, \f]
/// with \f$ f(0)=f'(0)=0 \f$ and \f$ f'(\infty) = 1 \f$
/// via Runge-Kutta and (scalar) Newton iteration.
///
/// Solving a BVP via local methods such as these is generally
/// not a great way of doing things. Actually, in this case
/// the BVP can formally be converted to an IVP and iteration
/// can be avoided altogether. Here we check the iteration plus
/// IVP approach by comparing it with the Cauchy problem formulation.

#include <IVP_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the Blasius equation by
    /// inheriting Equation base class.
    class Blasius_equation : public Equation<double>
    {
    public:

      /// Construct a 3rd order real system
      Blasius_equation( ) : Equation<double>( 3 ) {}

      /// We implement the equation as 3 first-order ODEs.
      void residual_fn( const DenseVector<double>& z, DenseVector<double>& f ) const
      {
        f[ 0 ] = z[ 1 ];
        f[ 1 ] = z[ 2 ];
        f[ 2 ] = -z[ 0 ] * z[ 2 ];
      }

      void mass( const DenseVector<double>&x, DenseMatrix<double> &m ) const
      {
        Utility::fill_identity(m);
      }
    };

    /// Define a residual function using the boundary
    /// conditions for the Blasius profile.
    class Blasius_residual : public Residual<double>
    {
    public:
      ODE_IVP<double> *ode;

      Blasius_residual() : Residual<double> ( 1 )
      {
        // instantiate the Blasius equation as our problem
        eqn = new Blasius_equation;
        // instantiate an ODE in [0,20] with a maximum of
        // 1000 integrator steps, using this problem equation
        ode = new ODE_IVP<double>( eqn, 0.0, 20.0, 1000 );
        // let's keep the history of the IVP for output in main
      }

      ~Blasius_residual()
      {
        // clean up the created objects
        delete ode;
        delete eqn;
      }

      /// implement a residual function.
      void residual_fn( const DenseVector<double>& unknown, DenseVector<double>& BC ) const
      {
        DenseVector<double> u( 3, 0.0 ); // 3 variables of Blasius problem
        u[ 0 ] = 0.0;          // f = 0
        u[ 1 ] = 0.0;          // f' = 0
        u[ 2 ] = unknown[ 0 ];      // f'' = unknown
        // shoot with the RKF method, tolerance of 1.e-9
        u = ode -> shoot45( u, 1.e-7, 0.01 );
        // return the residual of f'(infty) = 1 condition
        BC[ 0 ] = u[ 1 ] - 1.;
      }
    private:
      Blasius_equation *eqn;
    };
  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== BVP_Shoot: scalar Newton solution of Blasius ====\n";
  cout << "\n";

  Example::Blasius_residual problem;
  // instantiate a scalar Newton class for the one Blasius residual
  Newton<double> Blasius( &problem );

  DenseVector<double> stress( 1, 0.0 );
  stress[ 0 ] = 1.0;
  try
  {
    Blasius.iterate( stress );
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }

  // output the data
  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  TrackerFile my_file( "./DATA/Blasius.dat" );
  my_file.push_ptr( &( problem.ode -> get_mesh() ) );
  my_file.update();

  // There is a rescaling that removes any need to iterate
  // c = f'(infty) when stress = 1 at the plate
  const double c = 1.65519036023e0;
  const double answer = 1. / pow( c, 3. / 2. );

  const double tol( 1.e-6 );          // A pass/fail tolerance
  if ( abs( stress[0] - answer ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout.precision( 14 );
    cout << abs( stress[0] - answer ) << " > " << tol << "\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
}
