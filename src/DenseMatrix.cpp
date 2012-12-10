/// \file DenseMatrix.cpp
/// Implementation of a DENSE matrix as
/// an Vector of NVectors, inheriting from Matrix_base.

#include <complex>
#include <algorithm>

#include <Matrix_base.h>
#include <DenseVector.h>
#include <DenseMatrix.h>
#include <Exceptions.h>
#include <Functors.h>
#include <Utility.h>

namespace CppNoddy
{

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix( ) : Matrix_base<_Type>(), NR( 0 ), NC( 0 )
  {}

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix( const std::size_t& rows,
                                   const std::size_t& cols,
                                   const _Type& fill ) :
      Matrix_base<_Type>(),
      NR( rows ),
      NC( cols )
  {
    // make a row
    const DenseVector<_Type> row( cols, fill );
    // reserve the space
    MATRIX.reserve( rows );
    for ( std::size_t i = 0; i < rows; ++i )
    {
      // push require number of rows into the 'matrix'
      MATRIX.push_back( row );
    }
  }

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix( const std::size_t& rows,
                                   const std::size_t& cols,
                                   const _Type* p ) :
      Matrix_base<_Type>(),
      NR( rows ),
      NC( cols )
  {
    MATRIX.reserve( rows );
    for ( std::size_t i = 0; i < rows; ++i )
    {
      MATRIX.push_back( DenseVector<_Type>( cols, &p[ i * cols ] ) );
    }
  }

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix( const DenseMatrix<_Type>& source ) :
      Matrix_base<_Type>( source )
  {
    *this = source;
  }

  template <typename _Type>
  DenseMatrix<_Type>::~DenseMatrix( )
  {}

  template <typename _Type>
  inline DenseMatrix<_Type>& DenseMatrix<_Type>::operator=( const DenseMatrix<_Type>& source )
  {
    if ( this == &source )
      return * this;
    MATRIX = source.MATRIX;
    NR = source.NR;
    NC = source.NC;
    return *this;
  }

  template <typename _Type>
  std::size_t DenseMatrix<_Type>::nelts() const
  {
    return NR * NC;
  }

  template <typename _Type>
  void DenseMatrix<_Type>::add( const DenseMatrix<_Type>& B )
  {
#ifdef PARANOID
    // check number of columns at least match
    if ( ( B.nrows() != NR ) || ( B.ncols() != NC ) )
    {
      std::string problem( "The DenseMatrix.add has a geometry error.\n" );
      throw ExceptionGeom( problem, NR, NC, B.nrows(), B.ncols() );
    }
#endif
    std::transform( MATRIX.begin(), MATRIX.end(), B.MATRIX.begin(), MATRIX.begin(), std::plus< DenseVector<_Type> >() );
  }

  template <typename _Type>
  void DenseMatrix<_Type>::sub( const DenseMatrix<_Type>& B )
  {
#ifdef PARANOID
    // check number of columns at least match
    if ( ( B.nrows() != NR ) || ( B.ncols() != NC ) )
    {
      std::string problem( "The DenseMatrix.sub has a geometry error.\n" );
      throw ExceptionGeom( problem, NR, NC, B.nrows(), B.ncols() );
    }
#endif
    std::transform( MATRIX.begin(), MATRIX.end(), B.MATRIX.begin(), MATRIX.begin(), std::minus< DenseVector<_Type> >() );
  }

  template <typename _Type>
  void DenseMatrix<_Type>::scale( const _Type& mult )
  {
    std::transform( MATRIX.begin(), MATRIX.end(), MATRIX.begin(), scale_functor< DenseVector<_Type>, _Type >( mult ) );
  }

  template <typename _Type>
  void DenseMatrix<_Type>::transpose()
  {
    if ( nrows() == ncols() )
    {
      // square matrix needs no temp object
      // loop through upper half diagonal of the matrix
      for ( std::size_t i = 0; i < nrows(); ++i )
      {
        for ( std::size_t j = i + 1; j < ncols(); ++j )
        {
          // swap elements
          std::swap( MATRIX[ i ][ j ], MATRIX[ j ][ i ] );
        }
      }
    }
    else
    {
      std::vector< DenseVector<_Type> > temp;
      temp.resize( NC );
      for ( std::size_t row = 0; row < NC; ++row )
      {
        temp[ row ].resize( NR );
      }
      for ( std::size_t i = 0; i < nrows(); ++i )
      {
        for ( std::size_t j = 0; j < ncols(); ++j )
        {
          temp[ j ][ i ] = MATRIX[ i ][ j ];
        }
      }
      MATRIX = temp;
      std::swap( NR, NC );
    }
  }

  template <typename _Type>
  DenseVector<_Type> DenseMatrix<_Type>::multiply( const DenseVector<_Type>& X ) const
  {
#ifdef PARANOID
    // check number of columns at least match
    if ( X.size() != NC )
    {
      std::string problem( "The DenseMatrix.multiply has a geometry error.\n" );
      throw ExceptionGeom( problem, NR, NC, X.size(), 1 );
    }
#endif
    DenseVector<_Type> temp;
    temp.reserve( NR );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      temp.push_back( Utility::dot( MATRIX[ row ], X ) );
    }
    return temp;
  }

  template <typename _Type>
  DenseMatrix<_Type> DenseMatrix<_Type>::multiply( const DenseMatrix<_Type>& B ) const
  {
#ifdef PARANOID
    // check number of columns at least match
    if ( B.nrows() != NC )
    {
      std::string problem( "The DenseMatrix.multiply has a geometry error.\n" );
      throw ExceptionGeom( problem, NR, NC, B.nrows(), B.ncols() );
    }
#endif
    // temporary object for the result
    DenseMatrix<_Type> C( nrows(), B.ncols(), 0.0 );
    // loops thru the columns in the B matrix
    for ( std::size_t col_in_B = 0; col_in_B < B.ncols(); ++col_in_B )
    {
      // set the column in the result to be the matrix-vector
      // product of (*this).multiply( column in B )
      C.set_col( col_in_B , multiply( B.get_col( col_in_B ) ) );
    }
    return C;
  }


  template <typename _Type>
  typename std::vector<DenseVector<_Type> >::iterator
  DenseMatrix<_Type>::max_in_col( const std::size_t& col,
                                  row_iter row_min, row_iter row_max ) const
  {
    row_iter index( row_min );
    double maxelt( std::abs( *( row_min -> begin() ) ) );
    for ( row_iter row = row_min + 1; row != row_max ; ++row )
    {
      const double elt( std::abs( *( row -> begin() + col ) ) );
      if ( elt >= maxelt )
      {
        maxelt = elt;
        index = row;
      }
    }
    return index;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::one_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].one_norm() );
    }
    return max;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::two_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].two_norm() );
    }
    return max;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::inf_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].inf_norm() );
    }
    return max;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::frob_norm() const
  {
    double sum ( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      sum += MATRIX[ row ].two_norm();
    }
    return sum;
  }


  template <typename _Type>
  void DenseMatrix<_Type>::dump() const
  {
    std::cout << "DENSE mtx size = " << nrows() << " x " << ncols() << "; \n";
    std::cout.precision( 4 );
    std::cout << "- start matrix \n";
    for ( std::size_t i = 0; i < nrows(); ++i )
    {
      std::cout << " row " << i << " =  ";
      for ( std::size_t j = 0; j < ncols(); ++j )
      {
        std::cout << MATRIX[ i ][ j ] << ", ";
      }
      std::cout << "\n";
    }
    std::cout << "- end matrix \n";
  }

  template <typename _Type>
  void DenseMatrix<_Type>::set_col( const std::size_t& col, const DenseVector<_Type>& X )
  {
    for ( std::size_t row = 0; row < NR; ++row )
    {
      MATRIX[ row ][ col ] = X[ row ];
    }
  }

  template <typename _Type>
  DenseVector<_Type> DenseMatrix<_Type>::get_col( const std::size_t& col ) const
  {
    DenseVector<_Type> X( NR, 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      X[ row ] = MATRIX[ row ][ col ];
    }
    return X;
  }

  template <>
  void DenseMatrix<std::complex<double> >::matrix_to_vector( DenseVector<double> &p, const std::size_t &padding ) const
  {
    p.reserve( 2 * NR * NC );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      for ( std::size_t col = 0; col < NC; ++col )
      {
        p.push_back( MATRIX[ row ][ col ].real() );
        p.push_back( MATRIX[ row ][ col ].imag() );
      }
      for ( std::size_t col = 0; col < padding; ++col )
      {
        p.push_back( 0.0 );
        p.push_back( 0.0 );
      }
    }
  }


  template <>
  void DenseMatrix<double>::matrix_to_vector( DenseVector<double> &p, const std::size_t &padding ) const
  {
    p.reserve( NR * NC );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      for ( std::size_t col = 0; col < NC; ++col )
      {
        p.push_back( MATRIX[ row ][ col ] );
      }
      for ( std::size_t col = 0; col < padding; ++col )
      {
        p.push_back( 0.0 );
      }
    }
  }

  template <>
  DenseVector<double> DenseMatrix<double>::matrix_to_vector( const std::size_t &padding ) const
  {
    DenseVector<double> V;
    V.reserve( NR * NC );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      for ( std::size_t col = 0; col < NC; ++col )
      {
        V.push_back( MATRIX[ row ][ col ] );
      }
      for ( std::size_t col = 0; col < padding; ++col )
      {
        V.push_back( 0.0 );
      }
    }
    return V;
  }

  template <>
  DenseVector<double> DenseMatrix<std::complex<double> >::matrix_to_vector( const std::size_t &padding ) const
  {
    DenseVector<double> V;
    V.reserve( 2 * NR * NC );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      for ( std::size_t col = 0; col < NC; ++col )
      {
        V.push_back( MATRIX[ row ][ col ].real() );
        V.push_back( MATRIX[ row ][ col ].imag() );
      }
      for ( std::size_t col = 0; col < padding; ++col )
      {
        V.push_back( 0.0 );
        V.push_back( 0.0 ); 
      }
    }
    return V;
  }


  // the versions to be used are:
  template class DenseMatrix<double>
  ;
  template class DenseMatrix<std::complex<double> >
  ;

} // end namespace

