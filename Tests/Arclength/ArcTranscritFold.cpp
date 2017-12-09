/// \file ArcTranscritFold.cpp
/// \ingroup Tests
/// \ingroup Arclength
/// Solve the nonlinear scalar residual problem
/// \f[ R \equiv (x - 3)( (x-2)^2 + (p-4) ) = 0 \f]
/// by arc-length continuation from the starting solution
/// \f$ x=0, p=0 \f$. There is a limit point \f$ x=2, p =4 \f$ and a
/// transcritical bifurcation \f$ x =3, p=3 \f$ on the computed branch of
/// solutions, both of which should be flagged through an exception being
/// raised. The Example will fail if two bifurcations are not found.

#include <cassert>

#include <Newton_bundle.h>

namespace CppNoddy
{
  namespace Example
  {
    /// Define the residual for arc-length continuation
    class Arc_problem : public Residual<double>
    {
    public:
      double p;
      double eps;

      Arc_problem() : Residual<double>( 1 ) {}

      void residual_fn( const DenseVector<double> &z, DenseVector<double> &f ) const
      {
        f[0] = ( z[0] - 3 ) * ( ( z[0] - 2 ) * ( z[0] - 2 )
                                + ( p - 4 ) ) + eps;
      }
    };

  } // end Example namespace
} // end CppNoddy namespace


using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== ARC: Scalar continuation and bifn detection =====\n";
  cout << "\n";

  // Instantiate the problem
  Example::Arc_problem residual_problem;
  residual_problem.p = 0.0;
  residual_problem.eps = 0.0;

  const double tol = 1.e-10;
  // Scalar Newton iteration problem
  Newton<double> newton( &residual_problem, 8, tol );
  newton.set_monitor_det( true );
  // initial guess
  DenseVector<double> state( 1, 0.1 );
  // initialise a state for arc length cont.
  newton.init_arc( state, &residual_problem.p, 0.01, 0.5 );

  // output for the data
  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );
  TrackerFile my_file( "./DATA/Trans_Fold.dat" );
  my_file.push_ptr( &residual_problem.p, "parameter" );
  my_file.push_ptr( &state, "system state" );
  my_file.header();
  my_file.precision( 8 );

  unsigned n_bif( 0 );
  for ( int i = 0; i < 100; ++i )
  {
    try
    {
      try
      {
        newton.arclength_solve( state );
      }
      catch ( std::runtime_error )
      {
        cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
        return 1;
      }
    }
    catch ( ExceptionBifurcation )
    {
      cout << " Bifurcation detected near x = " << state[0] <<
           " p = " << residual_problem.p << "\n\n";
      ++n_bif;
    }
    my_file.update();
  }

  if ( n_bif != 2 )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    return 1;
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
    return 0;
  }
}
