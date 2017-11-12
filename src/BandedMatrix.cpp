/// \file BandedMatrix.cpp
/// Implementation of the matrix class that constructs a BANDED matrix as
/// an NVector of NVectors.

#include <string>
#include <memory>
#include <algorithm>

#include <DenseVector.h>
#include <Matrix_base.h>
#include <BandedMatrix.h>
#include <Exceptions.h>

namespace CppNoddy
{

  template <typename _Type>
  BandedMatrix<_Type>::BandedMatrix( const std::size_t& rows,
                                     const std::size_t& upper_offdiag_bands,
                                     const _Type& fill ) :
      Matrix_base<_Type>(),
      N( rows ),
      L( upper_offdiag_bands )
  {
    // we'll assume that the upper & lower offdiags are of equal size.
    // logical STORAGE is the main diagonal + upper offdiagonal + lower offdiagonal
    // = 2 * L + 1.
    // However, allowing for row interchanging during pivotting leads to additional
    // padding of L, to make a total of 3*L+1.
    STORAGE = DenseVector<_Type>( N * ( 3 * L + 1 ), 0.0 );
  }

  template <typename _Type>
  BandedMatrix<_Type>::BandedMatrix( const BandedMatrix<_Type>& source ) :
      Matrix_base<_Type>( source )
  {
    *this = source;
  }

  template <typename _Type>
  inline BandedMatrix<_Type>& BandedMatrix<_Type>::operator=( const BandedMatrix<_Type>& source )
  {
    if ( this == &source )
      return * this;
    STORAGE = source.STORAGE;
    N = source.N;
    L = source.L;
    return *this;
  }

  template <typename _Type>
  std::size_t BandedMatrix<_Type>::nelts() const
  {
    return STORAGE.size();
  }

  template <typename _Type>
  void BandedMatrix<_Type>::scale( const _Type& mult )
  {
    STORAGE *= mult;
  }


  template <typename _Type>
  DenseVector<_Type> BandedMatrix<_Type>::multiply( const DenseVector<_Type>& X ) const
  {
    std::string problem;
    problem = " The multiply method of the BandedMatrix class has \n";
    problem += " not been implemented. \n";
    throw ExceptionRuntime( problem );
    return STORAGE; // dummy return
  }


  template <typename _Type>
  void BandedMatrix<_Type>::transpose()
  {
    for ( std::size_t row = 0; row < N; ++row )
    {
      for ( std::size_t col = std::max( 0, int( row ) - int( L ) ); col < row; ++col )
      {
        // swap elements
        std::swap( operator()( row, col ), operator()( col, row ) );
      }
    }
  }


  template <typename _Type>
  void BandedMatrix<_Type>::row_swap( const std::size_t& row1, const std::size_t& row2 )
  {
    // WE MUST HAVE row2 > row1 & row1 must have no elements left of the diagonal!
    // only of any use to native Gaussian elimination solver
    /// \todo MAKE ME PRIVATE OR ELSE!
    for ( std::size_t col = row1; col <= std::min( row1 + 2 * L, N - 1 ); ++col )
    {
      std::swap( operator()( row1, col ), operator()( row2, col ) );
    }
  }


  template <typename _Type>
  double BandedMatrix<_Type>::one_norm() const
  {
    double max( 0.0 );
    for ( unsigned row = 0; row < N; ++row )
    {
      // copy the row into a temp vector
      DenseVector<_Type> temp( 3 * L + 1, &STORAGE[ row * ( 3 * L + 1 ) ] );
      max = std::max( max, temp.one_norm() );
    }
    return max;
  }


  template <typename _Type>
  double BandedMatrix<_Type>::two_norm() const
  {
    double max( 0.0 );
    for ( unsigned row = 0; row < N; ++row )
    {
      // copy the row into a temp vector
      DenseVector<_Type> temp( 3 * L + 1, &STORAGE[ row * ( 3 * L + 1 ) ] );
      max = std::max( max, temp.two_norm() );
    }
    return max;
  }


  template <typename _Type>
  double BandedMatrix<_Type>::inf_norm() const
  {
    double max( 0.0 );
    for ( unsigned row = 0; row < N; ++row )
    {
      // copy the row into a temp vector
      DenseVector<_Type> temp( 3 * L + 1, &STORAGE[ row * ( 3 * L + 1 ) ] );
      max = std::max( max, temp.inf_norm() );
    }
    return max;
  }


  template <typename _Type>
  double BandedMatrix<_Type>::frob_norm() const
  {
    double sum( 0.0 );
    for ( unsigned row = 0; row < N; ++row )
    {
      // copy the row into a temp vector
      DenseVector<_Type> temp( 3 * L + 1, &STORAGE[ row * ( 3 * L + 1 ) ] );
      sum += temp.two_norm();
    }
    return sum;
  }

  template <typename _Type>
  void BandedMatrix<_Type>::dump() const
  {
    std::cout << "BANDED mtx size = " << nrows() << " x " << ncols()
              << " with INPUT no. of bands = " << 2 * L + 1
              << " and total bands STORED (for pivotting) = " << 3 * L + 1 << "\n";
    std::cout.precision( 4 );
    std::cout << "- start matrix \n";
    for ( std::size_t i = 0; i < nrows(); ++i )
    {
      std::cout << " row " << i << " =  ";
      for ( std::size_t j = 0; j < nrows(); ++j )
      {
        if ( i == j )
          std::cout << "*";
        if ( ( ( int ) j < ( int ) i - ( int ) L ) || ( ( int ) j > ( int ) i + ( int ) ( 2*L ) ) )
        {
          std::cout << 0 << ", ";
        }
        else
        {
          std::cout << operator() ( i, j ) << ", ";
        }
      }
      std::cout << "\n";
    }
    std::cout << "- end matrix \n";
  }

  template <>
  double* BandedMatrix<double>::base()
  {
    return &( STORAGE[0] );
  }

  template <>
  double* BandedMatrix<std::complex<double> >::base()
  {
    //return &( STORAGE[0].real() );
    return &reinterpret_cast<double(&)[2]>( STORAGE[0] )[0];
  }

  // the templated versions we require are:
  template class BandedMatrix<double>
  ;
  template class BandedMatrix<std::complex<double> >
  ;

} // end namespace
