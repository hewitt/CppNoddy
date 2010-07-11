/// \file Utility.h
/// A spec for a collection of utility functions.

#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <ctime>
#include <numeric>
#include <algorithm>
#include <numeric>

#include <Types.h>
#include <Exceptions.h>
#include <Functors.h>
#include <Matrix_base.h>
#include <FortranBLAS.h>

namespace CppNoddy
{
  /// Some utility methods associated with CppNoddy containers.
  namespace Utility
  {

    // MISC UTILS

    /// Seed the random number generator using current clock
    void time_seed();

    //
    //
    // MATRIX/VECTOR FILLING UTILS
    //
    //

    /// Fill a specified entire row of a dense matrix
    /// \param A The dense matrix to be used
    /// \param row The row in A to be modified
    /// \param value The value to be written to the entire row
    template <typename _Type>
    void fill_row( DenseMatrix<_Type>& A, const std::size_t& row, const _Type& value )
    {
      for ( std::size_t col = 0; col < A.ncols(); ++col )
      {
        A( row, col ) = value;
      }
    }

    /// Fill a specified entire row of a banded matrix.
    /// \param A The banded matrix to be used
    /// \param row The row in A to be modified
    /// \param value The value to be written to the entire row
    template <typename _Type>
    void fill_row( BandedMatrix<_Type>& A, const std::size_t& row, const _Type& value )
    {
      std::size_t L = A.noffdiag();
      // here we only write to the whole STORED banded matrix
      // including those elements that are present purely for
      // partial pivotting.
      for ( std::size_t col = std::max( int( row ) - int( L ), 0 );
            col <= std::min( row + 2 * L, A.ncols() - 1 ); ++col )
      {
        A( row, col ) = value;
      }
    }

    /// Fill a diagonal band of a matrix
    /// \param A The matrix to be used
    /// \param offset The offset of the band from the main diagonal
    ///  e.g. 0 = main diagional, -1 = first sub-diagonal
    /// \param value The value to be written to the band elements
    template <typename _Type>
    void fill_band( Matrix_base<_Type>& A, const int& offset, const _Type& value )
    {
      for ( std::size_t row = 0; row < A.nrows(); ++row )
      {
        if ( ( row + offset < A.ncols() ) && ( row + offset >= 0 ) )
        {
          A( row, row + offset ) = value;
        }
      }
    }

    /// Set all elements of a specified BANDED matrix
    /// \param A The BANDED matrix to be filled
    /// \param value The value to be placed in each element of the banded matrix
    template <typename _Type>
    void fill( BandedMatrix<_Type>& A, const _Type& value )
    {
      for ( std::size_t row = 0; row < A.nrows(); ++row )
      {
        fill_row( A, row, value );
      }
    }

    /// Set all elements of a specified DENSE matrix
    /// \param A The DENSE matrix to be filled
    /// \param value The value to be placed in each element of the matrix
    template <typename _Type>
    void fill( DenseMatrix<_Type>& A, const _Type& value )
    {
      for ( std::size_t row = 0; row < A.nrows(); ++row )
      {
        fill_row( A, row, value );
      }
    }

    /// Set all elements of a DENSE vector
    /// \param X The DENSE vector to be filled
    /// \param value The value to be placed in each element of the vector
    template <typename _Type>
    void fill( DenseVector<_Type>& X, const _Type& value )
    {
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        X[ i ] = value;
      }
    }

    /// Fill the three main diagonals with given data
    /// Other bands are left unchanged.
    /// \param A The matrix to be filled
    /// \param L The lower diagonal data
    /// \param D The main diagonal data
    /// \param U The upper diagonal data
    template <typename _Type>
    void fill_tridiag( Matrix_base<_Type>& A, const _Type& L,
                       const _Type& D, const _Type& U )
    {
      // fill the lower, diagonal and upper bands
      fill_band( A, -1, L );
      fill_band( A, 0, D );
      fill_band( A, 1, U );
    }

    /// Fill the main diagonal with ones
    /// Other matrix elements are left unchanged
    /// \param A The matrix to be filled
    template <typename _Type>
    void fill_identity( Matrix_base<_Type>& A )
    {
      if ( A.nrows() != A.ncols() )
      {
        std::string problem;
        problem = " In Utilities::fill_identity you are trying to ";
        problem += " fill a non square matrix as an identity matrix. \n";
        throw ExceptionRuntime( problem );
      }
      fill_band( A, 0, ( _Type ) 1.0 );
    }

    /// Fill a DENSE vector with random entries
    /// \param V The DENSE vector to be filled
    void fill_random( DenseVector<double>& V );

    /// Fill a DENSE vector with COMPLEX entries
    /// \param V The complex DENSE vector to be filled
    void fill_random( DenseVector< D_complex >& V );

    /// Fill a SPARSE vector with a set number of real random
    /// entries.
    /// \param V The SPARSE vector to be filled
    /// \param num_of_elts The total number of elements required
    ///  to be filled
    void fill_random( SparseVector<double>& V, const unsigned& num_of_elts );

    /// Fill a SPARSE vector with a set number of complex random
    /// entries.
    /// \param V The SPARSE vector to be filled
    /// \param num_of_elts The total number of elements required
    ///  to be filled
    void fill_random( SparseVector<D_complex>& V, const unsigned& num_of_elts );

    /// Fill a DENSE matrix with random elements
    /// \param A The DENSE matrix to be filled
    void fill_random( DenseMatrix<double>& A );

    /// Fill a BANDED matrix with random data
    /// \param A The BANDED matrix to be filled
    void fill_random( BandedMatrix<double>& A );

    /// Return a DENSE vector with the nodal points of a uniform
    /// mesh distributed between the upper/lower bounds as specified
    /// \param lower The lower bound of the uniform nodal distribution
    /// \param upper The upper bound of the uniform nodal distribution
    /// \param N The number of nodal points
    DenseVector<double> uniform_node_vector( const double& lower, const double& upper, const std::size_t& N );

    /// Return a DENSE vector with the nodal points of a non-uniform
    /// mesh distributed between the upper/lower bounds as specified
    /// with more nodes clustered near lower or upper depending upon
    /// the differencee of the power from unity. When power=1 this should
    /// provide a uniform mesh.
    /// \param lower The lower bound of the uniform nodal distribution
    /// \param upper The upper bound of the uniform nodal distribution
    /// \param N The number of nodal points
    /// \param power A measure of the non-uniformity
    DenseVector<double> power_node_vector( const double& lower, const double& upper, const std::size_t& N, const double& power );

    /// Return a dense vector with two uniform distributions in two separate
    /// regions.
    /// \param lower The first node
    /// \param mid The node that defines the boundary between the uniform meshes
    /// \param upper The final node
    /// \param N1 The number of nodes in the first region
    /// \param N2 The number of nodes in the second region
    DenseVector<double> two_uniform_node_vector( const double& lower, const double& mid, const double& upper, const std::size_t& N1, const std::size_t& N2 );

    /// Return a dense vector of nodal positions with more nodes concentrated
    /// at the mid point of the range.
    /// \param lower The first nodal position.
    /// \param upper The final nodal position.
    /// \param N The number of nodes required.
    /// \param power A measure of the non-uniformity, power = 1 => uniform distribution
    DenseVector<double> mid_weighted_node_vector( const double& lower, const double& upper, const std::size_t& N, const double& power );

    //
    //
    // SOME typical ops
    //
    //

    /// Given a DenseMatrix<double> of a streamfunction, this will compute the velocities
    /// at the same nodal points using 2nd-order finite differencing assuming a
    /// uniform nodal point distribution in a Cartesian coordinate system.
    /// \param source A pointer to the dense matrix containing the stream function
    /// \param dx The x-step size.
    /// \param dy The y-step size.
    /// \param u The matrix to return the "u" velocities.
    /// \param v The matrix to return the "v" velocities.
    void vels_from_streamfn_Cartesian( const DenseMatrix<double>& source, const double& dx, const double& dy, DenseMatrix<double>& u, DenseMatrix<double>& v );

    /// Given a DenseMatrix<double> of a streamfunction, this will compute the velocities
    /// at the same nodal points using 2nd-order finite differencing assuming a
    /// uniform nodal point distribution in a cylindrical-polar coordinate system,
    /// where the left-most boundary is ASSUMED to be the axis of symmetry.
    /// \param source A pointer to the dense matrix containing the stream function
    /// \param dr The r-step size.
    /// \param dz The z-step size.
    /// \param u The matrix to return the "u" velocities.
    /// \param w The matrix to return the "w" velocities.
    void vels_from_streamfn_Stokes( const DenseMatrix<double>& source, const double& dr, const double& dz, DenseMatrix<double>& u, DenseMatrix<double>& w );

    //
    //
    // MATRIX OPERATIONS
    //
    //

    /// BLAS wrapper to do DOUBLE DENSE A_{MxK} * B_{KxN} = C_{MxN}
    /// Since this is a Fortran library, it assumes a column_major
    /// format, but CppNoddy uses row_major. To be consistent we'll
    /// simply do (B^T)_{NxK} * (A^T)_{KxM} = (C^T)_{NxM} instead.
    /// Note that inversion of the transpose of the result C^T
    /// is handled implicitly via the construction of C.
    /// \param A First dense double matrix to be multiplied
    /// \param B Second dense double matrix to be multiplied
    /// \return The result of the multiplication C=A*B
    DenseMatrix<double> multiply( DenseMatrix<double>& A, DenseMatrix<double>& B );

    //
    //
    // VECTOR OPERATIONS
    //
    //

    /// Templated dot product.
    /// \param X First dense vector
    /// \param Y Second dense vector
    /// \return The dot product
    template <typename _Type>
    _Type dot( const DenseVector<_Type>& X, const DenseVector<_Type>& Y )
    {
      if ( X.size() != Y.size() )
      {
        std::string problem;
        problem = "The Utilities::dot method has been called \n";
        problem += "with two unequal length vectors.";
        throw ExceptionGeom( problem, X.size(), Y.size() );
      }
      return inner_product( X.begin(), X.end(), Y.begin(), _Type( 0.0 ) );
    }

    template <typename _Type>
    int sgn( const _Type& a )
    {
      if ( a > ( _Type )0 )
      {
        return 1;
      }
      else
        if ( a < ( _Type )0 )
        {
          return -1;
        }
        else
        {
          return 0;
        }
    }

    //
    //
    // COMPLEX UTILS
    //
    //

    /// Return a double DENSE vector containing the real part
    /// of a complex DENSE vector
    /// \param X The complex vector to take the real part of
    DenseVector<double> real( const DenseVector<D_complex>& X );

    /// Return a double DENSE vector containing the imaginary part
    /// of a complex DENSE vector
    /// \param X The complex vector to take the imaginary part of
    DenseVector<double> imag( const DenseVector<D_complex>& X );

    //
    //
    // STRING TWEAKERY
    //
    //

    /// Return an integer value as a string - useful for file naming
    /// \param val The integer value to be stringified.
    std::string stringify( const int &val );

    /// Return a double value as a string - useful for file naming.
    /// \param val The double value to be stringified
    /// \param p Precision to be used in the output
    std::string stringify( const double &val, int p );

  }

} // end namespace

#endif // UTILITY_H
