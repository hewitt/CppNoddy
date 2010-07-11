/// \file BandedMatrix.h
/// A matrix class that constructs a BANDED matrix.
/// Storage is contiguous for use with external banded solvers.

#ifndef BANDEDMATRIX_H
#define BANDEDMATRIX_H

#include <DenseVector.h>
#include <Matrix_base.h>
#include <Exceptions.h>

namespace CppNoddy
{

  /// A matrix class that constructs a BANDED matrix.
  template <typename _Type>
  class BandedMatrix  : public Matrix_base<_Type>
  {

  public:

    typedef typename DenseVector<_Type>::elt_iter elt_iter;

    /// Empty constructor .. should we stop this?
    BandedMatrix()
    {}

    /// Noddy Banded Matrix constructor.
    /// \param rows The number of rows in the matrix.
    /// \param offdiag The maximum number of bands above OR below the diagonal.
    /// The total number stored will be 3 * offdiag + 1; i.e. assumes the
    /// number of offdiagonal elts is the same both above and below, and
    /// we must an extra 'offdiag' because of fill-in during pivotting.
    /// \param fill The initial value to be placed in each element
    ///   of the banded matrix.
    BandedMatrix( const std::size_t& rows, const std::size_t& offdiag, const _Type& fill );

    /// Copy constructor.
    /// \param source The source object to be copied
    BandedMatrix( const BandedMatrix& source );

    /// Assignment operator.
    /// \param source The source object for the assignment
    /// \return The newly assigned object
    BandedMatrix& operator=( const BandedMatrix& source );

    /// Access operator
    const _Type& operator() ( const std::size_t& row, const std::size_t& col ) const;
    /// Access operator
    _Type& operator() ( const std::size_t& row, const std::size_t& col );
    /// Access operator
    const _Type& get( const std::size_t& row, const std::size_t& col ) const;
    /// Access operator
    _Type& set( const std::size_t& row, const std::size_t& col );

    /// Get the number of rows
    /// \return The number of rows in the matrix
    std::size_t nrows() const;

    /// Get the number of columns
    /// \return The number of columns in the matrix
    std::size_t ncols() const;

    /// Get the number of elements
    /// \return The number of elements
    std::size_t nelts() const;

    /// Scale all entries in the matrix by a scalar
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

    /// Right multiply the matrix by a DENSE vector
    /// \param X The DENSE vector to multiply by
    /// \return A DENSE vector of length ncols() produced from
    ///  the multiplication
    DenseVector<_Type>  multiply( const DenseVector<_Type>& X ) const;

    /// Output the matrix contents to std::cout
    void dump() const;

    //
    // NON-INHERITED MEMBER FUNCTIONS
    //

    /// Assign a value the matrix so that it has the same
    /// geometry, but zero entries in all locations
    /// including those reserved for pivotting.
    /// \param elt The value to be assigned to all entries
    void assign( _Type elt )
    {
      STORAGE.assign( STORAGE.size(), elt );
    }

    /// Get the number of off-diagonal elements
    /// where the total INPUT band width is 2*noffdiag+1
    /// since the band structure is symmetric.
    /// \return The number of upper or lower offdiagonal bands;
    ///   NOTE: the number of offdiagonals STORED will typically
    ///  be three times this to allow for pivotting
    std::size_t noffdiag() const;

    /// Exchange rows in the matrix
    /// \param row1 First row to be swapped
    /// \param row2 Second row to be swapped
    void row_swap( const std::size_t& row1, const std::size_t& row2 );

    /// Allow direct access to the vector STORAGE.
    /// Dangerous, but used for passing to LAPACK.
    double* base();

    elt_iter get_elt_iter( std::size_t row, std::size_t col )
    {
      return STORAGE.begin() + L * ( 3 * col + 2 ) + row;
    }

    //private:

    /// A contiguous vector
    DenseVector<_Type> STORAGE;
    /// The number of rows/cols in the matrix.
    std::size_t N;
    /// Max number of (INPUT) bands above OR below the main diagonal
    std::size_t L;

  }
  ; // end class


  // INLINED METHODS ARE BELOW

  template <typename _Type >
  inline const _Type& BandedMatrix<_Type>::operator() ( const std::size_t& row, const std::size_t& col ) const
  {
#ifdef PARANOID
    // if outside the NxN matrix
    if ( ( row >= N ) || ( row < 0 ) || ( col >= N ) || ( col < 0 ) )
    {
      std::string problem;
      problem = " The const operator() of BandedMatrix has been called \n";
      problem += " with a (row, col) index that is outside \n";
      problem += " the square NxN matrix.\n";
      throw ExceptionGeom( problem, N, N, row, col );
    }
    // check if the subscripts are out of the band
    if ( !( ( col + L >= row ) && ( col <= 2*L + row ) ) )
    {
      std::string problem;
      problem = " The const operator() of BandedMatrix has been called \n";
      problem += " with a (row, col) index that is outside \n";
      problem += " the band structure. Bandwidth and offset from\n";
      problem += " the diagonal information follows.\n";
      throw ExceptionGeom( problem, L, col - row );
    }
#endif
    // MOSTLY WE PUSH BANDED MATRICES TO LAPACK - so we'll keep column major
    // internal storage ... a la Fortran
    // return STORAGE[ col * ( 3 * L + 1 ) + ( row - col ) + 2 * L ];
    // or equiv.
    return STORAGE[ L * ( 3 * col + 2 ) + row ];
  }

  template <typename _Type>
  inline _Type& BandedMatrix<_Type>::operator() ( const std::size_t& row, const std::size_t& col )
  {
#ifdef PARANOID
    // if outside the NxN matrix
    if ( ( row >= N ) || ( row < 0 ) || ( col >= N ) || ( col < 0 ) )
    {
      std::string problem;
      problem = " The operator() of BandedMatrix has been called \n";
      problem += " with a (row, col) index that is outside \n";
      problem += " the square NxN matrix.\n";
      throw ExceptionRange( problem, N, N, row, col );
    }
    // check if the subscripts are out of the band
    if ( !( ( col + L >= row ) && ( col <= 2*L + row ) ) )
    {
      std::string problem;
      problem = " The operator() of BandedMatrix has been called \n";
      problem += " with a (row, col) index that is outside \n";
      problem += " the band structure.\n";
      std::cout << " L = " << L << "\n";
      std::cout << " row = " << row << "\n";
      std::cout << " col = " << col << "\n";
      throw ExceptionRange( problem, N, L, row, col );
    }
#endif
    // MOSTLY WE PUSH BANDED MATRICES TO LAPACK - so we'll keep column major
    // internal storage ... a la Fortran
    // return STORAGE[ col * ( 3 * L + 1 ) + ( row - col ) + 2 * L ];
    // or equiv.
    return STORAGE[ L * ( 3 * col + 2 ) + row ];
  }

  template <typename _Type>
  inline const _Type& BandedMatrix<_Type>::get
  ( const std::size_t& row, const std::size_t& col ) const
  {
    return operator() ( row, col );
  }

  template <typename _Type>
  inline _Type& BandedMatrix<_Type>::set
  ( const std::size_t& row, const std::size_t& col )
  {
    return operator() ( row, col );
  }

  template <typename _Type>
  inline std::size_t BandedMatrix<_Type>::nrows() const
  {
    return N;
  }

  template <typename _Type>
  inline std::size_t BandedMatrix<_Type>::ncols() const
  {
    return N;
  }

  template <typename _Type>
  inline std::size_t BandedMatrix<_Type>::noffdiag() const
  {
    return L;
  }



} // end namespace

#endif
