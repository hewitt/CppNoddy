/// \file FnQuadrature.cpp
/// Implementation of the real quadrature class. Note that
/// the object uses a function pointer for the integrand.

#include <string>

#include <Types.h>
#include <FnQuadrature.h>
#include <Exceptions.h>
#include <Utility.h>

namespace CppNoddy
{

  // ctor
  FnQuadrature::FnQuadrature( fn_ptr ptr_to_fn,
                              const double& x1,
                              const double& x2,
                              const unsigned& num_of_regions ) :
      A( x1 ),
      B( x2 ),
      p_FN( ptr_to_fn )
  {
    NODES = Utility::uniform_node_vector( x1, x2, num_of_regions );
  }

  // ctor
  FnQuadrature::FnQuadrature( fn_ptr ptr_to_fn,
                              const double& x1,
                              const double& x2,
                              const DenseVector<double>& nodes ) :
      A( x1 ),
      B( x2 ),
      p_FN( ptr_to_fn )
  {
    NODES = nodes;
  }

  void FnQuadrature::set_subintervals( const unsigned& n )
  {
    NODES = Utility::uniform_node_vector( A, B, n );
  }

  double FnQuadrature::Gauss( const int& n )
  {
    DenseVector<double> t;    // evaluation points in [-1,1]
    DenseVector<double> w;    // respective weight values
    switch ( n )
    {
    case 1:
      t.push_back( 0.0 );
      w.push_back( 2.0 );
      break;
    case 2:
      t.push_back( -1. / sqrt( 3. ) );
      w.push_back( 1.0 );
      t.push_back( 1. / sqrt( 3. ) );
      w.push_back( 1.0 );
      break;
    case 3:
      t.push_back( 0.0 );
      w.push_back( 8. / 9. );
      t.push_back( -0.774596669241483 );
      w.push_back( 5. / 9. );
      t.push_back( 0.774596669241483 );
      w.push_back( 5. / 9. );
      break;
    default:
      std::string problem;
      problem = " The Quadrature.Gauss method is trying to apply \n";
      problem += " a Gauss rule with more points than 3. \n";
      problem += " Currenlty only n=1,2,3 are supported! \n";
      throw ExceptionRuntime( problem );
    }

    double sum = 0.0;
    double c = 0.5 * ( B + A );
    double m = 0.5 * ( B - A );
    double f, x;
    t.scale( m );

    for ( unsigned i = 0; i < t.size(); ++i )
    {
      x = c + t[ i ];
      p_FN( x , f );
      sum += w[ i ] * f;
    }

    sum = sum * m;

    return sum;
  }


  double FnQuadrature::sub_Gauss( const int& n )
  {
    double x_left , x_right;
    double sum = 0.0;
    unsigned num = NODES.size();
    //    double h = ( B - A ) / num;

    for ( unsigned i = 0; i < num - 1; ++i )
    {
      //x_left = A + i * h;
      x_left = NODES[ i ];
      //x_right = x_left + h;
      x_right = NODES[ i + 1 ];
      FnQuadrature I( p_FN, x_left , x_right , num );
      sum += I.Gauss( n );
    }

    return sum;
  }

  double FnQuadrature::trapezium()
  {
    double sum = 0.0;
    unsigned num = NODES.size();
    double h = ( B - A ) / num;
    double x_left , x_right;
    double f_left , f_right;

    x_left = A;
    p_FN( x_left , f_left );

    for ( unsigned i = 0; i < num; ++i )
    {
      x_right = A + ( i + 1 ) * h;
      p_FN( x_right , f_right );
      sum += ( f_left + f_right );
      f_left = f_right;
    }
    sum *= h / 2.;
    return sum;
  }


}

