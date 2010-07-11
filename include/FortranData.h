/// \file FortranData.h

#ifndef FORTRANDATA_H
#define FORTRANDATA_H

#include <Types.h>
#include <Uncopyable.h>
#include <Exceptions.h>

namespace CppNoddy
{
  /// A little (legacy) utility class for passing CppNoddy containers to
  /// Fortran library routines.
  class FortranData : private Uncopyable
  {

  public:

    /// Make an empty FortranData object that is not constructed
    /// from a CppNoddy container
    explicit FortranData( std::size_t size )
    {
      p_DATA = new DenseVector<double>( size, 0.0 );
    }

    /// Make a FortranData object from a double dense matrix (with potential)
    /// for padding the data to avoid cache contention. The transpose flag
    /// defaults to true in order to return data in column-major format.
    /// \param matrix A double dense matrix to convert to Lapack/Fortran format
    /// \param transpose Transpose the matrix to be in column major format
    /// \param padding The number of additional elements to pad the matrix
    FortranData( DenseMatrix<double>& matrix, bool transpose = true, int padding = 0 )
    {
      p_DATA = new DenseVector<double>;
      if ( transpose )
      {
        matrix.transpose();
      }
      matrix.matrix_to_vector( *p_DATA, padding );
    }

    /// Make a FortranData object from a complex dense matrix (with potential)
    /// for padding the data to avoid cache contention. The transpose flag
    /// defaults to true in order to return data in column-major format.
    /// \param matrix A complex dense matrix to convert to Lapack/Fortran format
    /// \param transpose Transpose the matrix to be in column major format
    /// \param padding The number of additional elements to pad the matrix
    FortranData( DenseMatrix<D_complex>& matrix, bool transpose = true, int padding = 0 )
    {
      p_DATA = new DenseVector<double>;
      if ( transpose )
      {
        matrix.transpose();
      }
      matrix.matrix_to_vector( *p_DATA, padding );
    }

    FortranData( BandedMatrix<double>& matrix, bool transpose = true )
    {
      std::string problem;
      problem = "The default storage format for a banded matrix is now\n";
      problem += "a contiguous column-major vector. You should not need to \n";
      problem += "the FortranData method anymore.\n";
      throw ExceptionRuntime( problem );
    }

    FortranData( BandedMatrix<D_complex>& matrix, bool transpose = true )
    {
      std::string problem;
      problem = "The default storage format for a banded matrix is now\n";
      problem += "a contiguous column-major vector. You should not need to \n";
      problem += "the FortranData method anymore.\n";
      throw ExceptionRuntime( problem );
    }

    FortranData( DenseVector<D_complex>& vector )
    {
      std::string problem;
      problem = "This method is not required. Just use the base address of \n";
      problem += "the first element and treat it as a vector of doubles.\n";
      throw ExceptionRuntime( problem );
    }

    ~FortranData()
    {
      delete p_DATA;
    }

    /// Get the pointer to the first element.
    /// \return The pointer to the base address.
    double* base()
    {
      return & ( p_DATA -> operator[] ( 0 ) );
    }

    /// Get the reference to the first element.
    /// \return The the i-th element in the vector
    double& operator[] ( const std::size_t& i )
    {
      return p_DATA -> operator[] ( i );
    }

    /// Find the number of stored elements.
    /// \return The number of elements.
    std::size_t size() const
    {
      return p_DATA -> size();
    }

    /// Convert the data to a double dense format.
    /// \param rows The number of required rows
    /// \param cols The number of required columns
    /// \return The DenseMatrix<double> of the data
    DenseMatrix<double> to_dense_matrix( std::size_t rows, std::size_t cols )
    {
      DenseMatrix<double> matrix( rows, cols, &( p_DATA->operator[] ( 0 ) ) );
      return matrix;
    }

    /// Dump the data to standard out.
    void dump()
    {
      p_DATA -> dump();
    }

  private:

    // pointer to the base address of the data
    DenseVector<double>* p_DATA;
  };
}

#endif /*FORTRANDATA_H*/
