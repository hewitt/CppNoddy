/// \file SparseMatrix.h
/// A matrix class that constructs a SPARSE matrix as
/// an STL Vector of SparseVectors, inheriting from Matrix_base.

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <SparseVector.h>
#include <DenseVector.h>
#include <Exceptions.h>
#include <Sequential_Matrix_base.h>

#if defined(PETSC_D) || defined(PETSC_Z)
  #include "petsc.h"
#endif

namespace CppNoddy {

  template <typename _SystemType>
  class SparseLinearSystem;

  /// A matrix class that constructs a SPARSE matrix as a
  /// row major std::vector of SparseVectors.
  template <typename _Type>
  class SparseMatrix : public Sequential_Matrix_base<_Type> {
    typedef typename std::map< std::size_t, _Type >::const_iterator citer;
    typedef typename std::map< std::size_t, _Type >::iterator iter;

   public:

    /// Construct with a set number of rows
    /// \param rows The number of rows in the matrix
    /// \param cols The number of columns in the matrix
    SparseMatrix(const std::size_t& rows, const std::size_t& cols);

    /// Construct from a row permutation of another sparse matrix
    /// \param source_rows Defines the permutation, row i of this matrix
    /// is row source_rows[i] of the source
    SparseMatrix(const SparseMatrix<_Type>& source, const std::vector<std::size_t>& source_rows);

    /// Copy constructor.
    /// \param source The source object to be copied
    SparseMatrix(const SparseMatrix& source);

    /// Assignment operator.
    /// \param source The source object for the assignment
    /// \return The newly assigned object
    SparseMatrix& operator=(const SparseMatrix& source);

    /// Default d-tor
    ~SparseMatrix() {}

    /// Blank the contents of this matrix
    void blank() {
      m_matrix.clear();
      m_matrix.reserve(m_nr);
      SparseVector<_Type> sparse_row(m_nc);
      for(std::size_t i = 0; i < m_nr; ++i) {
        m_matrix.push_back(sparse_row);
      }
    }

    /// Access operator
    const _Type& operator()(const std::size_t& row, const std::size_t& col) const;
    /// Access operator
    _Type& operator()(const std::size_t& row, const std::size_t& col);
    /// Access operator
    const _Type& get(const std::size_t& row, const std::size_t& col) const;
    /// Access operator
    _Type& set(const std::size_t& row, const std::size_t& col);

    
    /// Get a row of the matrix
    /// \param row The index of the row to be get
    /// \return A sparse vector of values in the row
    SparseVector<_Type> get_row(const std::size_t& row) const {
      return m_matrix[row];
    }

    /// Set a row of the matrix
    /// \param row The index of the row to be set
    /// \param A sparse vector of values to go into the row
    void set_row(const std::size_t& row, const SparseVector<_Type>& row_vector) {
      m_matrix[row] = row_vector;
    }

    /// Get the number of rows in the matrix
    /// \return The number of rows
    std::size_t nrows() const;

    /// Get the number of columns in the matrix
    /// \return The number of columns
    std::size_t ncols() const;

    /// Get the number of (non-zero) elements in the matrix
    /// \return The number of (non-zero) elements
    std::size_t nelts() const;

    /// Scale the matrix by a scalar
    /// \param mult The scalar multiplier
    void scale(const _Type& mult);

    /// Transpose the matrix in place
    // void transpose();
    
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
    /// \return A DENSE vector of the result of the multiplication
    DenseVector<_Type> multiply(const DenseVector<_Type>& X) const;
    

    /// Output the contents of the matrix to std::cout
    void dump() const;

    /// A simple method for dumping the matrix to a file
    /// \param filename The filename to write the data to (will overwrite)
    /// \param precision Precision of the output strings
    void dump(std::string filename, int precision = 10) const {
      std::ofstream dump;
      dump.open(filename.c_str());
      dump.precision(precision);
      dump.setf(std::ios::showpoint);
      dump.setf(std::ios::showpos);
      dump.setf(std::ios::scientific);
      for(std::size_t row = 0; row < m_nr; ++row) {
        for(citer pos = m_matrix[row].begin(); pos != m_matrix[row].end(); ++pos) {
          dump << row << " " << pos -> first << " " << pos -> second << "\n";
        }
      }
      dump << "\n";
      dump.close();
    }

    //
    // NON-INHERITED MEMBER FUNCTIONS
    //

    /// Operator overloading for ROW access
    /// \param row The row to access
    /// \return The DENSE vector of the row data
    SparseVector<_Type>& operator[](const std::size_t& row);

    /// Operator overloading for ROW access
    /// \param row The row to access
    /// \return The DENSE vector of the row data
    const SparseVector<_Type>& operator[](const std::size_t& row) const;


#if defined (PETSC_D) || defined (PETSC_Z)
    // PETSc is compiled separately for double or complex ... but only linked against
    // once in the build process, so it's either/or.

    /// Takes the contents of the SparseMatrix and converts it to a
    /// standard compressed format for a specified row.
    /// \param row_number The row index to extract the data for
    /// \param storage A contiguous vector of the non-zero elts
    ///  (has to be allocated and big enough)
    /// \param cols The column indices of each entry in the storage vector
    ///  (has to be allocated and big enough)
    void get_row_petsc(PetscInt row_number, PetscScalar* storage, PetscInt* cols);
    
    /// Extracts the number of non-zero elements in each row and returns them as a
    /// PetscInt array of length m_nr. The array should be allocated on entry.
    /// \param The array to fill with the number of non-zero elts in each row
    ///  (has to be allocated and big enough)
    void nelts_all_rows(PetscInt* row_nnz) {
      for(std::size_t i = 0; i < m_nr; ++i) {
        row_nnz[i] = m_matrix[i].nelts();
      }
    }
#endif

    /// The number of non-zero elements in a specified row
    /// \param row The row index to return the number of non-zero elts for
    std::size_t nelts_in_row(int row) {
      return m_matrix[row].nelts();
    }

    /// Find the maximum entry in a column -- used in the native solver.
    /// \param col The column to search through
    /// \param row_min The start row for the search
    /// \param row_max The end row for the search (NOT INCLUSIVE)
    std::size_t max_in_col(const std::size_t& col, const std::size_t& row_min,
                           const std::size_t& row_max) const;

    /// Swap two rows in the matrix -- used in the native solver.
    /// \param row1 The first row to be exchanged
    /// \param row2 The second row to be exchanged
    void row_swap(const std::size_t& row1, const std::size_t& row2);

   private:

    /// An STL vector of SparseVectors.
    std::vector< SparseVector<_Type> > m_matrix;
    /// The max number of rows in the matrix.
    std::size_t m_nr;
    /// The max number of columns in the matrix.
    std::size_t m_nc;

    template <typename _SystemType>
    friend class SparseLinearSystem;

  }
  ; // END CLASS


  // INLINED METHODS FOLLOW

  template <typename _Type>
  inline const _Type& SparseMatrix<_Type>::operator()(const std::size_t& row, const std::size_t& col) const {
    return m_matrix[ row ].get(col);
  }

  template <typename _Type>
  inline _Type& SparseMatrix<_Type>::operator()(const std::size_t& row, const std::size_t& col) {
    return m_matrix[ row ][ col ];
  }

  template <typename _Type>
  inline const _Type& SparseMatrix<_Type>::get(const std::size_t& row, const std::size_t& col) const {
    return this -> operator()(row, col);
  }

  template <typename _Type>
  inline _Type& SparseMatrix<_Type>::set(const std::size_t& row, const std::size_t& col) {
    return this -> operator()(row, col);
  }

  template <typename _Type>
  inline SparseVector<_Type>& SparseMatrix<_Type>::operator[](const std::size_t& row) {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The SparseMatrix.get_row has a range error.\n");
      throw ExceptionRange(problem, m_nr, row);
    }
#endif
    return m_matrix[ row ];
  }


  template <typename _Type>
  inline const SparseVector<_Type>& SparseMatrix<_Type>::operator[](const std::size_t& row) const {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The SparseMatrix.get_row has a range error.\n");
      throw ExceptionRange(problem, m_nr, row);
    }
#endif
    return m_matrix[ row ];
  }

  template <typename _Type>
  inline std::size_t SparseMatrix<_Type>::nrows() const {
    return m_nr;
  }

  template <typename _Type>
  inline std::size_t SparseMatrix<_Type>::ncols() const {
    return m_nc;
  }

  template <typename _Type>
  inline void SparseMatrix<_Type>::row_swap(const std::size_t& row1, const std::size_t& row2) {
    // actually do the swap, rather than keep a row permutation vector.
    // std::swap<SparseVector<_Type> > ( matrix[ row1 ], matrix[ row2 ] );
    m_matrix[ row1 ].swap(m_matrix[ row2 ]);
  }


} // end namespace


#endif
