/// \file Equation_with_double_mass.cpp
/// Implementation for an equations class that can be inherited from
/// to allow instantiation of ODE objects using the resulting class.

#include <Equation_with_double_mass.h>
#include <Residual_with_coords.h>
#include <Utility.h>

namespace CppNoddy
{

  template <typename _Type>
  Equation_with_double_mass<_Type >::Equation_with_double_mass( const unsigned& order ) :
      Residual_with_coords<_Type>( order, 3 )
  {
    // initialise the container for the mass matrices
    MASS1_AT_LAST_STATE = DenseMatrix<_Type>( order, order, 0.0 );
    MASS2_AT_LAST_STATE = DenseMatrix<_Type>( order, order, 0.0 );
  }

  template <typename _Type>
  Equation_with_double_mass<_Type >::~Equation_with_double_mass()
  {
    // timer reporting is done via the Equation (base) class
  }

  template <typename _Type>
  void Equation_with_double_mass<_Type >::update( const DenseVector<_Type> &state )
  {
    // call the base class's update method
    Residual_with_coords<_Type>::update( state );
    // this has to go after the base class update - otherwise we'll be
    // resuming the timer twice in a row
#ifdef TIME
    this -> T_UPDATER.start();
#endif
    // now deal with the two mass matrices separately
    mass1( state, MASS1_AT_LAST_STATE );
    mass2( state, MASS2_AT_LAST_STATE );
#ifdef TIME
    this -> T_UPDATER.stop();
#endif
  }

  template <typename _Type>
  void Equation_with_double_mass<_Type >::get_jacobian_of_mass1_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const
  {
    // we dont need state in the default implementation as its already been set by the update method
    //
    // copy some items for FD computation of Jacobian of mass matrix
    DenseVector<_Type> copy_of_state( this -> LAST_STATE );
    DenseMatrix<_Type> copy_of_mass( MASS1_AT_LAST_STATE );
    std::vector< DenseMatrix<_Type> > jacmass;
    // update the Jacobian of the mass matrix
    for ( std::size_t i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      copy_of_state[ i ] += this -> DELTA;
      mass1( copy_of_state, copy_of_mass );
      copy_of_state[ i ] -= this -> DELTA;
      copy_of_mass.sub( MASS1_AT_LAST_STATE );
      copy_of_mass.scale( 1. / this -> DELTA );
      // the 3D object that represents the Jacobian of the mass matrix
      jacmass.push_back( copy_of_mass );
    }
    // evaluate the jacabian of mass contribution
    for ( unsigned i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      for ( unsigned j = 0; j < this -> ORDER_OF_SYSTEM; ++j )
      {
        h( i, j ) = Utility::dot( jacmass[ j ][ i ], vec );
      }
    }
  }

  template <typename _Type>
  void Equation_with_double_mass<_Type>::get_jacobian_of_mass2_mult_vector( const DenseVector<_Type> &state, const DenseVector<_Type> &vec, DenseMatrix<_Type> &h ) const
  {
    // we dont need state in the default implementation as its already been set by the update method
    //
    // copy some items for FD computation of Jacobian of mass matrix
    DenseVector<_Type> copy_of_state( this -> LAST_STATE );
    DenseMatrix<_Type> copy_of_mass( MASS2_AT_LAST_STATE );
    std::vector< DenseMatrix<_Type> > jacmass;
    // update the Jacobian of the mass matrix
    for ( std::size_t i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      copy_of_state[ i ] += this -> DELTA;
      mass2( copy_of_state, copy_of_mass );
      copy_of_state[ i ] -= this -> DELTA;
      copy_of_mass.sub( MASS2_AT_LAST_STATE );
      copy_of_mass.scale( 1. / this -> DELTA );
      // the 3D object that represents the Jacobian of the mass matrix
      jacmass.push_back( copy_of_mass );
    }
    // evaluate the jacabian of mass contribution
    for ( unsigned i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
    {
      for ( unsigned j = 0; j < this -> ORDER_OF_SYSTEM; ++j )
      {
        h( i, j ) = Utility::dot( jacmass[ j ][ i ], vec );
      }
    }
  }

  // the required templated versions are:
  template class Equation_with_double_mass<double>
  ;
  template class Equation_with_double_mass<std::complex<double> >
  ;

} // end namespace
