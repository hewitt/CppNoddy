/// \file Quad.cpp
/// \ingroup Examples
/// \ingroup Generic
/// Compute the integral
/// \f[ \int_0^{20} \cos(x) \exp \left ( -\frac{x}{4} \right ) \mbox{d}x \f] with varying
/// schemes, then compare the result to the exact value
/// \f[ \frac{4 e^{-5}}{17} \left ( e^5 - \cos (20) + 4\sin (20) \right ) \f]

#include <Generic_bundle.h>
#include <FnQuadrature.h>

namespace CppNoddy
{
  namespace Example
  {
    bool failed = false;

    /// The function that defines the integrand.
    void Qfn( const double &x, double &f )
    {
      f = std::cos( x ) * std::exp( -x / 4. );
    }

    /// The test will fail if it requires more than a set number
    /// of subintervals in the integration to converge to the
    /// exact value to within the specified tol.
    void test( unsigned n )
    {
      if ( n > 16777216 )
      {
        failed = true;
      }
    }

  } // end Example namespace
} // end CppNoddy namespace

using namespace CppNoddy;
using namespace std;

int main()
{

  cout << "\n";
  cout << "=== Quadrature: integral of cos(x)exp(-x/4) =========\n";
  cout << "\n";

  /// Define I as being a Quadrature object with limits 0,20
  /// and a default of 1 sub-interval.
  FnQuadrature I( Example::Qfn, 0.0, 20.0 , 1 );

  double result;
  unsigned n = 1;
  const double exact = ( 4. / 17. ) * ( exp( 5.0 ) - cos( 20.0 ) + 4 * sin( 20.0 ) ) * exp( -5.0 );
  const double tol = 1.e-7;
  cout << " Tolerance used : " << tol << "\n";
  cout << "\n";
  cout << " Trapezium summation \n";
  cout.precision( 12 );

  do
  {
    I.set_subintervals( n );
    result = I.trapezium() - exact;
#ifdef DEBUG
    cout << "n = " << n << "  Integral error = " << abs( result ) << "\n";
#endif

    n *= 2;
    Example::test( n );
  }
  while ( ( abs( result ) > tol ) && !Example::failed );

  cout << "    required " << n / 2 << " sub-intervals. \n";
  cout << "\n";
  cout << " Sub_Gauss with 1 Guass point \n";

  n = 1;
  do
  {
    I.set_subintervals( n );
    result = I.sub_Gauss( 1 ) - exact;
#ifdef DEBUG
    cout << "n = " << n << "  Integral error = " << abs( result ) << "\n";
#endif

    n *= 2;
    Example::test( n );
  }
  while ( ( abs( result ) > tol ) && !Example::failed );

  cout << "    required " << n / 2 << " sub-intervals. \n";
  cout << "\n";
  cout << " Sub_Gauss with 2 Guass points \n";

  n = 1;
  do
  {
    I.set_subintervals( n );
    result = I.sub_Gauss( 2 ) - exact;
#ifdef DEBUG
    cout << "n = " << n << "  Integral error = " << abs( result ) << "\n";
#endif

    n *= 2;
    Example::test( n );
  }
  while ( ( abs( result ) > tol ) && !Example::failed );

  cout << "    required " << n / 2 << " sub-intervals. \n";
  cout << "\n";
  cout << " Sub_Gauss with 3 Guass points \n";

  n = 1;
  do
  {
    I.set_subintervals( n );
    result = I.sub_Gauss( 3 ) - exact;
#ifdef DEBUG
    cout << "n = " << n << "  Integral error = " << abs( result ) << "\n";
#endif

    n *= 2;
    Example::test( n );
  }
  while ( ( abs( result ) > tol ) && !Example::failed );

  cout << "    required " << n / 2 << " sub-intervals. \n";

  if ( Example::failed )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }
}
