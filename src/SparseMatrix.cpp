/// \file SparseMatrix.cpp
/// Implementation of a SPARSE matrix as
/// an STL Vector of SparseVectors, inheriting from Matrix_base.

#include <string>
#include <complex>

#include <SparseMatrix.h>
#include <Matrix_base.h>
#include <Exceptions.h>

namespace CppNoddy
{

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix( const std::size_t& rows, const std::size_t& cols )
      : Matrix_base<_Type>(), NR( rows ), NC( cols )
  {
    MATRIX.reserve( NR );
    SparseVector<_Type> sparse_row( NC );
    for ( std::size_t i = 0; i < NR; ++i )
    {
      MATRIX.push_back( sparse_row );
    }
  }

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix( const SparseMatrix<_Type>& source ) :
      Matrix_base<_Type>( source )
  {
    *this = source;
  }

  template <typename _Type>
  inline SparseMatrix<_Type>& SparseMatrix<_Type>::operator=( const SparseMatrix<_Type>& source )
  {
    if ( this == &source )
      return * this;
    MATRIX = source.MATRIX;
    NR = source.NR;
    NC = source.NC;
    return *this;
  }

  template <typename _Type>
  std::size_t SparseMatrix<_Type>::nelts() const
  {
    std::size_t num_of_elts( 0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      num_of_elts += MATRIX[ row ].nelts();
    }
    return num_of_elts;
  }

  template <typename _Type>
  std::size_t SparseMatrix<_Type>::max_in_col( const std::size_t& col,
      const std::size_t& row_min, const std::size_t& row_max ) const
  {
    double maxelt( 0.0 );
    // return outside of the array as default
    std::size_t index = nrows();
    for ( std::size_t row = row_min; row < row_max ; ++row )
    {
      // only bother looking up entries of rows with a first
      // element in a column less than the one we're looking at
      if ( MATRIX[ row ].begin() -> first <= col )
      {
        const double elt( std::abs( MATRIX[ row ].get( col ) ) );
        if ( elt >= maxelt )
        {
          maxelt = elt;
          index = row;
        }
      }
    }
    return index;
  }

  template <typename _Type>
  void SparseMatrix<_Type>::scale( const _Type& mult )
  {
    for ( std::size_t row = 0; row < NR; ++row )
    {
      MATRIX[ row ] *= mult;
    }
  }

  template <typename _Type>
  void SparseMatrix<_Type>::transpose()
  {
    throw ExceptionRuntime( "SparseMatrix.transpose has not been implemented" );
  }

  template <typename _Type>
  double SparseMatrix<_Type>::one_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].one_norm() );
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::two_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].two_norm() );
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::inf_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].inf_norm() );
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::frob_norm() const
  {
    double sum ( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      sum += MATRIX[ row ].two_norm();
    }
    return sum;
  }

  template <typename _Type>
  DenseVector<_Type> SparseMatrix<_Type>::multiply( const DenseVector<_Type>& X ) const
  {
    throw ExceptionRuntime( "SparseMatrix.multiply has not been implemented" );
  }

  template <typename _Type>
  void SparseMatrix<_Type>::dump() const
  {
    std::cout << "SPARSE mtx size = " << NR << " x  sparse \n";
    std::cout.precision( 4 );
    std::cout << "- start matrix \n";
    for ( std::size_t row = 0; row < NR; ++row )
    {
      std::cout << " row " << row << " :  ";
      MATRIX[ row ].dump();
      std::cout << " \n";
    }
    std::cout << "- end matrix \n";
  }

  // the versions to be used are:
  template class SparseMatrix<double>
  ;
  template class SparseMatrix<std::complex<double> >
  ;

} // end namespace

