/// \file Equation_1matrix.cpp
/// Implementation for an equations class that can be inherited from
/// to allow instantiation of ODE objects using the resulting class.

#include <Equation_1matrix.h>
#include <Residual_with_coords.h>
#include <Utility.h>

namespace CppNoddy
{

  template <typename _Type, typename _Xtype>
  Equation_1matrix<_Type, _Xtype>::Equation_1matrix( const unsigned& order ) :
      Residual_with_coords<_Type, _Xtype>( order, 1 )
  {
    MATRIX0_AT_LAST_STATE = DenseMatrix<_Type>( this -> ORDER_OF_SYSTEM, this -> ORDER_OF_SYSTEM, 0.0 );
  }

  template <typename _Type, typename _Xtype>
  Equation_1matrix<_Type, _Xtype>::~Equation_1matrix()
  {
    // timer reporting is done via the Equation class this inherits from
  }

  template <typename _Type, typename _Xtype>
  void Equation_1matrix<_Type, _Xtype>::get_jacobian_of_matrix0_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const
  {
    // we dont need state in the default implementation as its already been set by the update method. You do need it for the user
    // to overload this method with an explicit analytical version however.
    //
    // copy some items for FD computation of Jacobian of matrix
    DenseVector<_Type> copy_of_state( this -> LAST_STATE );
    DenseMatrix<_Type> copy_of_matrix( MATRIX0_AT_LAST_STATE );
    std::vector< DenseMatrix<_Type> > jacmatrix;
    // update the Jacobian of the mass matrix
    for ( std::size_t i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      copy_of_state[ i ] += this -> DELTA;
      matrix0( copy_of_state, copy_of_matrix );
      copy_of_state[ i ] -= this -> DELTA;
      copy_of_matrix.sub( MATRIX0_AT_LAST_STATE );
      copy_of_matrix.scale( 1. / this -> DELTA );
      // the 3D object that represents the Jacobian of the mass matrix
      jacmatrix.push_back( copy_of_matrix );
    }
    // evaluate the jacabian of mass contribution
    for ( unsigned i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      for ( unsigned j = 0; j < this -> ORDER_OF_SYSTEM; ++j )
      {
        h( i, j ) = Utility::dot( jacmatrix[ j ][ i ], vec );
      }
    }
  }

  template <typename _Type, typename _Xtype>
  void Equation_1matrix<_Type, _Xtype>::update( const DenseVector<_Type> &state )
  {
    // use the base class update to set LAST_STATE, FN_AT_LAST_STATE and JAC_AT_LAST_STATE
    Residual_with_coords<_Type, _Xtype>::update( state );
    // this has TIME ifdef to go after the base class update - otherwise we'll be
    // resuming the timer twice in a row
#ifdef TIME
    this -> T_UPDATER.start();
#endif
    // now we deal with the mass matrix separately to set MATRIX0_AT_LAST_STATE
    // and JAC_OF_MATRIX
    matrix0( state, MATRIX0_AT_LAST_STATE );
#ifdef TIME
    this -> T_UPDATER.stop();
#endif
  }

  // the required templated versions are:
  template class Equation_1matrix<double>
  ;
  template class Equation_1matrix<std::complex<double> >
  ;
  template class Equation_1matrix<std::complex<double>, std::complex<double> >
  ;

} // end namespace
