/// \file SparseMatrix.h
/// A matrix class that constructs a SPARSE matrix as
/// an STL Vector of SparseVectors, inheriting from Matrix_base.

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>

#include <SparseVector.h>
#include <DenseVector.h>
#include <Matrix_base.h>
#include <Exceptions.h>


namespace CppNoddy
{

  template <typename _SystemType>
  class SparseLinearSystem;

  /// A matrix class that constructs a SPARSE matrix as a
  /// row major std::vector of SparseVectors.
  template <typename _Type>
  class SparseMatrix : public Matrix_base<_Type>
  {
    typedef typename std::map< std::size_t, _Type >::const_iterator citer;
    typedef typename std::map< std::size_t, _Type >::iterator iter;

  public:

    /// Construct with a set number of rows
    /// \param rows The number of rows in the matrix
    /// \param cols The number of columns in the matrix
    SparseMatrix( const std::size_t& rows, const std::size_t& cols );

    /// Copy constructor.
    /// \param source The source object to be copied
    SparseMatrix( const SparseMatrix& source );

    /// Assignment operator.
    /// \param source The source object for the assignment
    /// \return The newly assigned object
    SparseMatrix& operator=( const SparseMatrix& source );

    ~SparseMatrix() {}

    /// Access operator
    const _Type& operator() ( const std::size_t& row, const std::size_t& col ) const;
    /// Access operator
    _Type& operator() ( const std::size_t& row, const std::size_t& col );
    /// Access operator
    const _Type& get( const std::size_t& row, const std::size_t& col ) const;
    /// Access operator
    _Type& set( const std::size_t& row, const std::size_t& col );

    /// Get the number of rows in the matrix
    /// \return The number of rows
    std::size_t nrows() const;

    /// Get the number of columns in the matrix
    /// \return The number of columns
    std::size_t ncols() const;

    /// Get the number of elements in the matrix
    /// \return The number of elements
    std::size_t nelts() const;

    /// Scale the matrix by a scalar
    /// \param mult The scalar multiplier
    void scale( const _Type& mult );

    /// Transpose the matrix
    void transpose();

    /// \return The maximum one_norm of all rows
    double one_norm() const;

    /// \return The maximum two_norm of all rows
    double two_norm() const;

    /// \return The maximum inf_norm of all rows
    double inf_norm() const;

    /// \return The sum of the two_norm of all rows
    double frob_norm() const;

    /// Right-multiply by a DENSE vector
    /// \param X The DENSE vector to be multiplied by
    /// \return A DENSE matrix of the result of the multiplication
    DenseVector<_Type> multiply( const DenseVector<_Type>& X ) const;

    /// Output the contents of the matrix to std::cout
    void dump() const;

    //
    // NON-INHERITED MEMBER FUNCTIONS
    //

    /// Operator overloading for ROW access
    /// \param row The row to access
    /// \return The DENSE vector of the row data
    SparseVector<_Type>& operator[] ( const std::size_t& row );

    /// Operator overloading for ROW access
    /// \param row The row to access
    /// \return The DENSE vector of the row data
    const SparseVector<_Type>& operator[] ( const std::size_t& row ) const;


    void get_row_compressed( _Type* storage, int* cols, int* rows )
    {
      // iterator to the maps that are used in SparseVector
      // this is bad form as it exposes the internals of the SparseVector storage
      citer pos;
      //std::size_t last_col( NC + 1 ); // no column should have NC + 1, so this is a dummy start value
      std::size_t i( 0 ); // where we are in the storage vector
      //
      for ( std::size_t row = 0; row < NR; ++row )
      {
        //std::cout << row << " \n";
        // flag to indicate that we're on a new coloumn
        bool new_row( true );
        pos = MATRIX[ row ].begin();
        do
        {
          _Type elt( pos -> second );
          int col( pos -> first );
          storage[ i ] = elt;
          cols[ i ] = col;
          ++pos;
          if ( new_row )
          {
            rows[ row ] = i;
            new_row = false;
          }
          ++i;
        }
        while ( pos != MATRIX[ row ].end() );
      }
      // last entry points to end + 1
      rows[ NR ] = nelts();
      //for ( std::size_t i = 0; i < nelts(); ++i )
      //{
      //std::cout << storage[ i ] << " ";
      //}
      //std::cout << "\n";
      //for ( std::size_t i = 0; i < nelts(); ++i )
      //{
      //std::cout << cols[ i ] << " ";
      //}
      //std::cout << "\n";
      //for ( std::size_t i = 0; i < NR + 1; ++i )
      //{
      //std::cout << rows[ i ] << " ";
      //}
      //std::cout << "\n";
    }

    /// Take the contents of the SparseMatrix and convert it to a
    /// standard compressed column format that can be fed directly
    /// into the SuperLU library to solve.
    /// \param storage A contiguous vector of the non-zero elts
    /// \param rows The row indices of each entry in the storage vector
    /// \param cols The indices of the elements in the storage vector
    ///   that begin a new column
    void get_col_compressed( _Type* storage, int* rows, int* cols )
    {
      // iterator to the maps that are used in SparseVector
      // this is bad form as it exposes the internals of the SparseVector storage
      citer pos;
      //std::size_t last_col( NC + 1 ); // no column should have NC + 1, so this is a dummy start value
      std::size_t i( 0 ); // where we are in the storage vector
      //
      for ( std::size_t col = 0; col < NC; ++col )
      {
        std::cout << col << " \n";
        // flag to indicate that we're on a new coloumn
        bool new_col( true );
        for ( std::size_t row = 0; row < NR; ++row )
        {
          // look for this element in the sparse matrix
          pos = MATRIX[ row ].find( col );
          // if found
          if ( pos != MATRIX[ row ].end() )
          {
            // row/col data exists
            storage[ i ] = pos -> second;
            rows[ i ] = row;
            if ( new_col )
            {
              // column starts at index i
              cols[ col ] = i;
              new_col = false;
            }
            ++i;
          }
        }
      }
      // last entry points to end + 1
      cols[ NC ] = nelts();
    }

  private:

    /// Find the maximum entry in a column
    /// \param col The column to search through
    /// \param row_min The start row for the search
    /// \param row_max The end row for the search (NOT INCLUSIVE)
    std::size_t max_in_col( const std::size_t& col, const std::size_t& row_min, const std::size_t& row_max ) const;

    /// Swap two rows in the matrix
    /// \param row1 The first row to be exchanged
    /// \param row2 The second row to be exchanged
    void row_swap( const std::size_t& row1, const std::size_t& row2 );

    /// An vector of SparseVectors.
    std::vector< SparseVector<_Type> > MATRIX;
    /// The max number of rows in the matrix.
    std::size_t NR;
    /// The max number of columns in the matrix.
    std::size_t NC;

    template <typename _SystemType>
    friend class SparseLinearSystem;

  }
  ; // END CLASS


  // INLINED METHODS FOLLOW

  template <typename _Type>
  inline const _Type& SparseMatrix<_Type>::operator() ( const std::size_t& row, const std::size_t& col ) const
  {
    return MATRIX[ row ].get( col );
  }

  template <typename _Type>
  inline _Type& SparseMatrix<_Type>::operator() ( const std::size_t& row, const std::size_t& col )
  {
    return MATRIX[ row ][ col ];
  }

  template <typename _Type>
  inline const _Type& SparseMatrix<_Type>::get( const std::size_t& row, const std::size_t& col ) const
  {
    return this -> operator() ( row, col );
  }

  template <typename _Type>
  inline _Type& SparseMatrix<_Type>::set( const std::size_t& row, const std::size_t& col )
  {
    return this -> operator() ( row, col );
  }

  template <typename _Type>
  inline SparseVector<_Type>& SparseMatrix<_Type>::operator[]( const std::size_t& row )
  {
#ifdef PARANOID
    if ( ( row > NR ) || ( row < 0 ) )
    {
      std::string problem( "The SparseMatrix.get_row has a range error.\n" );
      throw ExceptionRange( problem, NR, row );
    }
#endif
    return MATRIX[ row ];
  }


  template <typename _Type>
  inline const SparseVector<_Type>& SparseMatrix<_Type>::operator[]( const std::size_t& row ) const
  {
#ifdef PARANOID
    if ( ( row > NR ) || ( row < 0 ) )
    {
      std::string problem( "The SparseMatrix.get_row has a range error.\n" );
      throw ExceptionRange( problem, NR, row );
    }
#endif
    return MATRIX[ row ];
  }

  template <typename _Type>
  inline std::size_t SparseMatrix<_Type>::nrows() const
  {
    return NR;
  }

  template <typename _Type>
  inline std::size_t SparseMatrix<_Type>::ncols() const
  {
    return NC;
  }

  template <typename _Type>
  inline void SparseMatrix<_Type>::row_swap( const std::size_t& row1, const std::size_t& row2 )
  {
    // actually do the swap, rather than keep a row permutation vector.
    // std::swap<SparseVector<_Type> > ( matrix[ row1 ], matrix[ row2 ] );
    MATRIX[ row1 ].swap( MATRIX[ row2 ] );
  }


} // end namespace


#endif
