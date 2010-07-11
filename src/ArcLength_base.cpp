/// \file ArcLength_base.cpp
/// An implementation of the arclength class. This defines
/// all the usual get/set methods for arclength parameters and the
/// additional residual used when augmenting the system to allow
/// for iteration on the arclength.

#include <Uncopyable.h>
#include <ArcLength_base.h>

namespace CppNoddy
{

  template <typename _Type>
  double& ArcLength_base<_Type>::ds()
  {
    return DS;
  }

  template <typename _Type>
  double& ArcLength_base<_Type>::arcstep_multiplier()
  {
    return ARCSTEP_MULTIPLIER;
  }

  template <typename _Type>
  bool& ArcLength_base<_Type>::rescale_theta()
  {
    return RESCALE_THETA;
  }

  template <typename _Type>
  double& ArcLength_base<_Type>::theta()
  {
    return THETA;
  }

  template <typename _Type>
  double& ArcLength_base<_Type>::desired_arc_proportion()
  {
    return DESIRED_ARC_PROPORTION;
  }

  template <typename _Type>
  void ArcLength_base<_Type>::solve( DenseVector<_Type> &x )
  {}

  template <typename _Type>
  void ArcLength_base<_Type>::update( const DenseVector<_Type>& x )
  {
    if ( RESCALE_THETA )
    {
      update_theta( x );
    }
    X_DERIV_S = ( x - LAST_X ) / DS;
    PARAM_DERIV_S = ( *p_PARAM - LAST_PARAM ) / DS;
    LAST_X = x;
    LAST_PARAM = *p_PARAM;
  }

  template <typename _Type>
  void ArcLength_base<_Type>::init_arc( DenseVector<_Type> x,
                                        _Type* param,
                                        const double& length,
                                        const double& max_length )
  {
    // set the pointers to the parameter & state variable
    p_PARAM = param;
    // compute the solution at this parameter value in usual way
    solve( x );
    // store the converged solution at this point
    LAST_X = x;
    LAST_PARAM = *p_PARAM;
    //
    // we now have one solution and can get ready to arc-length continue
    //
    DS = length;
    // put the arc-length entirely into a parameter for step 1
    *p_PARAM += DS;
    // recompute the state at this parameter value
    solve( x );
    // update the derivatives of state & parameter variables
    update( x );
    INITIALISED = true;
    MAX_DS = max_length;
  }



  template <typename _Type>
  double ArcLength_base<_Type>::arclength_residual( const DenseVector<_Type>& x ) const
  {
    return THETA * ( x - LAST_X ).two_norm() / x.size()
           + ( 1.0 - THETA ) * std::pow( std::abs( *p_PARAM - LAST_PARAM ), 2 )
           - DS * DS;
  }

  template <typename _Type>
  void ArcLength_base<_Type>::update_theta( const DenseVector<_Type>& x )
  {
    if ( RESCALE_THETA )
    {
      double Delta_p2 = std::pow( std::abs( *p_PARAM - LAST_PARAM ), 2 );
      double Delta_x2 = ( x - LAST_X ).two_norm() / x.size();
      THETA = Delta_p2 * ( DESIRED_ARC_PROPORTION - 1.0 )
              / ( Delta_p2 * ( DESIRED_ARC_PROPORTION - 1.0 )
                  - DESIRED_ARC_PROPORTION * Delta_x2 );
    }
  }

  // the templated versions we require are:
  template class ArcLength_base<double>
  ;
  template class ArcLength_base<std::complex<double> >
  ;
}

