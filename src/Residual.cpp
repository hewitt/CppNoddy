/// \file Residual.cpp
/// Implementation of the (double/complex) VECTOR residual class. To be
/// inherited by anything that is to be passed to the Newton object.

#include <Residual.h>

namespace CppNoddy
{
  template <typename _Type>
  Residual<_Type>::Residual( const unsigned& order ) : DELTA( 1.e-8 )
  {
    // set the order of the vector system and initialise the member data
    ORDER_OF_SYSTEM = order;
    NUMBER_OF_VARS = order;
    LAST_STATE = DenseVector<_Type>( NUMBER_OF_VARS, 0.0 );
    FN_AT_LAST_STATE = DenseVector<_Type>( ORDER_OF_SYSTEM, 0.0 );
    JAC_AT_LAST_STATE = DenseMatrix<_Type>( ORDER_OF_SYSTEM, NUMBER_OF_VARS, 0.0 );
#ifdef TIME
    // timer
    T_UPDATER = Timer( "Updating of the residual object:" );
#endif
  }

  template <typename _Type>
  Residual<_Type>::Residual( const unsigned& order, const unsigned& nstate ) : DELTA( 1.e-8 )
  {
    // stet the order of the vector system and initialise the member data
    ORDER_OF_SYSTEM = order;
    NUMBER_OF_VARS = nstate;
    LAST_STATE = DenseVector<_Type>( NUMBER_OF_VARS, 0.0 );
    FN_AT_LAST_STATE = DenseVector<_Type>( ORDER_OF_SYSTEM, 0.0 );
    JAC_AT_LAST_STATE = DenseMatrix<_Type>( ORDER_OF_SYSTEM, NUMBER_OF_VARS, 0.0 );
#ifdef TIME
    // timer
    T_UPDATER = Timer( "Updating of the residual object:" );
#endif
  }

  template <typename _Type>
  Residual<_Type>::~Residual()
  {
#ifdef TIME
    std::cout << "\n";
    T_UPDATER.stop();
    T_UPDATER.print();
#endif
  }

  template <typename _Type>
  void Residual<_Type>::jacobian( const DenseVector<_Type>& state, DenseMatrix<_Type>& jac ) const
  {
    DenseVector<_Type> new_state( state );
    // evaluation of the function
    DenseVector<_Type> f_at_new_state( ORDER_OF_SYSTEM, 0.0 );
    // default is to FD the Jacobian
    for ( std::size_t i = 0; i < NUMBER_OF_VARS; ++i )
    {
      new_state[ i ] += DELTA;
      residual_fn( new_state, f_at_new_state );
      new_state[ i ] -= DELTA;
      jac.set_col( i, ( f_at_new_state - FN_AT_LAST_STATE ) / DELTA );
    }
  }

  // the templated versions we require are:
  template class Residual<double>
  ;
  template class Residual<std::complex<double> >
  ;

}   // end namespace
