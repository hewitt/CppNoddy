/// \file DenseMatrix.h
/// A matrix class that constructs a DENSE matrix as
/// an STL Vector of DenseVectors.

#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include <vector>
#include <Sequential_Matrix_base.h>
#include <DenseVector.h>

namespace CppNoddy {

  /// A matrix class that constructs a DENSE matrix as a
  /// row major std::vector of DenseVectors. This is generally a nice
  /// way to implement the matrix, but the data elements [i][j] are
  /// not necessarily in contiguous row_major format in memory
  /// because there is a system dependent padding at the end
  /// of each row vector. Thus, in general, &[i][j+1] - &[i][j] !=
  /// &[i][0] - &[i-1][Nc] .. ie. the step may be larger from the
  /// end of row i-1 to the start of row i. Typically, the data must
  /// be copied into contiguous memory for use in external libraries
  /// that take base pointers & assume uniform steps between elements.
  template <typename _Type>
  class DenseMatrix : public Sequential_Matrix_base<_Type> {
   public:

    /// Typedef iterator types
    typedef typename std::vector<DenseVector<_Type> >::iterator row_iter;
    typedef typename std::vector<DenseVector<_Type> >::const_iterator
    row_citer;
    typedef typename std::vector<DenseVector<_Type> >::reverse_iterator
    row_riter;
    typedef typename std::vector<DenseVector<_Type> >::const_reverse_iterator
    row_criter;
    typedef typename DenseVector<_Type>::elt_iter elt_iter;
    typedef typename DenseVector<_Type>::elt_citer elt_citer;
    typedef typename DenseVector<_Type>::elt_riter elt_riter;
    typedef typename DenseVector<_Type>::elt_criter elt_criter;

    /// Allow empty construction
    DenseMatrix();

    /// Noddy Matrix constructor with a specified fill data.
    /// \param rows The number of rows in the matrix.
    /// \param cols The number of columns in the matrix.
    /// \param fill The entry to be placed in all elements.
    DenseMatrix(const std::size_t & rows,
                const std::size_t & cols,
                const _Type & fill);

    /// Construct a Noddy Matrix from a contiguous set of data.
    /// This will be nasty if you pass the wrong pointer, but is
    /// useful in interfacing with external libraries.
    /// This assumes the contiguous data is in row_major format.
    /// \param rows The number of rows in the matrix.
    /// \param cols The number of columns in the matrix.
    /// \param p    A pointer to the start of the data.
    DenseMatrix(const std::size_t& rows,
                const std::size_t& cols,
                const _Type* p);

    /// Construct a dense matrix from its banded counterpart.
    /// \param source The banded matrix to be used in the construction.
    DenseMatrix(const BandedMatrix<_Type>& source) {
      // get the number of off diagonal elts
      int l = source.noffdiag();
      // banded matrix class is always square
      int n = source.nrows();
      for(int row = 0; row < n; ++row) {
        DenseVector<_Type> vecrow(n, 0.0);
        for(int col = std::max(row - l, 0);
            col < std::min(n, row + l + 1);
            ++col) {
          vecrow[ col ] = source(row, col);
        }
        m_matrix.push_back(vecrow);
      }
      m_nr = m_nc = n;
    }


    /// Copy constructor.
    /// \param source The source object to be copied
    DenseMatrix(const DenseMatrix& source);

    /// Assignment operator.
    /// \param source The source object for the assignment
    /// \return The newly assigned object
    DenseMatrix& operator=(const DenseMatrix& source);

    ~DenseMatrix();

    /// Access operator
    const _Type& operator()(const std::size_t& row, const std::size_t& col) const;
    /// Access operator
    _Type& operator()(const std::size_t& row, const std::size_t& col);
    /// Access operator
    const _Type& get(const std::size_t& row, const std::size_t& col) const;
    /// Access operator
    _Type& set(const std::size_t& row, const std::size_t& col);

    /// SIMD operator sugar
    DenseMatrix<_Type> operator*(const double& m) const {
      DenseMatrix<_Type> temp(*this);
      temp.scale(m);
      return temp;
    }
    DenseMatrix<_Type> operator+(const DenseMatrix<_Type>& A) const {
      DenseMatrix<_Type> temp(*this);
      for(std::size_t row = 0; row < m_nr; ++row) {
        for(std::size_t col = 0; col < m_nc; ++col) {
          temp(row, col) += A(row, col);
        }
      }
      return temp;
    }
    DenseMatrix<_Type> operator-(const DenseMatrix<_Type>& A) const {
      DenseMatrix<_Type> temp(*this);
      for(std::size_t row = 0; row < m_nr; ++row) {
        for(std::size_t col = 0; col < m_nc; ++col) {
          temp(row, col) -= A(row, col);
        }
      }
      return temp;
    }

    /// \return The number of rows
    std::size_t nrows() const;
    /// \return The number of columns
    std::size_t ncols() const;
    /// \return The number of elements
    std::size_t nelts() const;

    /// Scale all matrix elements by a scalar
    /// \param mult The scalar multiplier
    void scale(const _Type& mult);

    /// Transpose the matrix
    void transpose();

    /// Return the maximum one_norm of all rows
    double one_norm() const;

    /// Rreturn the maximum two_norm of all rows
    double two_norm() const;

    /// Return the maximum inf_norm of all rows
    double inf_norm() const;

    /// Return the sum of the two_norm of all rows
    double frob_norm() const;

    
    /// Right multiply the matrix by a DENSE vector
    /// \param x The DENSE vector to be multiplied
    /// \return The DENSE vector of the multiplication
    DenseVector<_Type> multiply(const DenseVector<_Type>& x) const;
    
    /// Output the matrix to std::cout
    void dump() const;

    //
    // NON-INHERITED MEMBER FUNCTIONS
    //

    /// Assign a value to the matrix but keep the same
    /// geometry.
    /// \param elt The value to be assigned to all entries
    void assign(_Type elt) {
      m_matrix.assign(m_matrix.size(), DenseVector<_Type>(m_nc, elt));
    }

    /// Pass through of iterator calls.
    /// \return Iterator to the begin row
    row_iter begin() {
      return m_matrix.begin();
    }

    /// Pass through of iterator calls.
    /// \return Iterator to the end row
    row_iter end() {
      return m_matrix.end();
    }

    /// Pass through of iterator calls.
    /// \return Reverse iterator to the begin row
    row_riter rbegin() {
      return m_matrix.rbegin();
    }

    /// Pass through of iterator calls.
    /// \return Reverse iterator to the end row
    row_riter rend() {
      return m_matrix.rend();
    }

    /// Pass through of iterator calls.
    /// \return Const iterator to the begin row
    row_citer begin() const {
      return m_matrix.begin();
    }

    /// Pass through of iterator calls.
    /// \return Const iterator to the end row
    row_citer end() const {
      return m_matrix.end();
    }

    /// Pass through of iterator calls.
    /// \return Const reverse iterator to the begin row
    row_criter rbegin() const {
      return m_matrix.rbegin();
    }

    /// Pass through of iterator calls.
    /// \return Const reverse iterator to the end row
    row_criter rend() const {
      return m_matrix.rend();
    }

    /// Operator overloading for ROW access
    /// \param row The row to access
    /// \return The DENSE vector of the row data
    DenseVector<_Type>& operator[](const std::size_t& row);

    /// Operator overloading for ROW access
    /// \param row The row to access
    /// \return The DENSE vector of the row data
    const DenseVector<_Type>& operator[](const std::size_t& row) const;

    /// Add a DENSE matrix to this object
    /// \param b The DENSE matrix to add to 'this'
    void add(const DenseMatrix<_Type>& b);

    /// Subtract a DENSE matrix from this object
    /// \param b The DENSE matrix to subtract from 'this'
    void sub(const DenseMatrix<_Type>& b);

    /// Right-multiply by a DENSE matrix and return the result
    /// \param b The matrix to right-multiply by
    /// \return The matrix result of the multiplication
    DenseMatrix<_Type> multiply(const DenseMatrix<_Type>& b) const;

    /// Set a column of the matrix
    /// \param col The column to be set
    /// \param x The DENSE vector of column information
    void set_col(const std::size_t& col, const DenseVector<_Type>& x);

    /// Get a column of the matrix
    /// \param col The column to get
    /// \return A DENSE vector of the column of data
    DenseVector<_Type> get_col(const std::size_t& col) const;

    /// Conversion to contiguous data in row major format
    /// Inefficient ... the void method is preferred
    /// \param padding An integer padding hack
    /// \return A DENSE vector containing the matrix
    DenseVector<double> matrix_to_vector(const std::size_t &padding = 0) const;

    /// Conversion to contiguous data in row major format
    /// \param padding An integer padding hack
    /// \param p A DENSE vector containing the matrix
    void matrix_to_vector(DenseVector<double> &p, const std::size_t &padding = 0) const;

    /// Find the maximum abs value in a column
    /// \param col The column offset to search through
    /// \param row_min Iterator to the begin row
    /// \param row_max Iterator to final row (NOT INCLUSIVE)
    /// \return An iterator to the row that contains the maximum value
    row_iter max_in_col(const std::size_t& col, row_iter row_min, row_iter row_max) const;

   private:

    /// An NVector of NVectors.
    std::vector< DenseVector<_Type> > m_matrix;
    /// The number of rows in the matrix.
    std::size_t m_nr;
    /// The number of columns in the matrix.
    std::size_t m_nc;

  }
  ; // END CLASS


  // INLINE METHODS BELOW


  template <typename _Type>
  inline const _Type& DenseMatrix<_Type>::operator()(const std::size_t& row,
      const std::size_t& col) const {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The DenseMatrix.() has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
    if((col > m_nc) || (col < 0)) {
      std::string problem("The DenseMatrix.() has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
#endif
    return m_matrix[ row ][ col ];
  }

  template <typename _Type>
  inline _Type& DenseMatrix<_Type>::operator()(const std::size_t& row,
      const std::size_t& col) {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The DenseMatrix.() has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
    if((col > m_nc) || (col < 0)) {
      std::string problem("The DenseMatrix.() has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
#endif
    return m_matrix[ row ][ col ];
  }

  template <typename _Type>
  inline const _Type& DenseMatrix<_Type>::get
  (const std::size_t& row, const std::size_t& col) const {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The DenseMatrix.get has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }

    if((col > m_nc) || (col < 0)) {
      std::string problem("The DenseMatrix.get has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
#endif
    return m_matrix[ row ][ col ];
  }

  template <typename _Type>
  inline _Type& DenseMatrix<_Type>::set
  (const std::size_t& row,
   const std::size_t& col) {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The DenseMatrix.set has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
    if((col > m_nc) || (col < 0)) {
      std::string problem("The DenseMatrix.set has a range error.\n");
      throw ExceptionRange(problem, m_nr, row, m_nc, col);
    }
#endif
    return m_matrix[ row ][ col ];
  }

  template <typename _Type>
  inline DenseVector<_Type>& DenseMatrix<_Type>::operator [](const std::size_t& row) {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The DenseMatrix.get_row has a range error.\n");
      throw ExceptionRange(problem, m_nr, row);
    }
#endif
    return m_matrix[ row ];
  }

  template <typename _Type>
  inline const DenseVector<_Type>& DenseMatrix<_Type>::operator [](const std::size_t& row) const {
#ifdef PARANOID
    if((row > m_nr) || (row < 0)) {
      std::string problem("The DenseMatrix.get_row has a range error.\n");
      throw ExceptionRange(problem, m_nr, row);
    }
#endif
    return m_matrix[ row ];
  }

  template <typename _Type>
  inline std::size_t DenseMatrix<_Type>::nrows() const {
    return m_nr;
  }

  template <typename _Type>
  inline std::size_t DenseMatrix<_Type>::ncols() const {
    return m_nc;
  }


} // end namespace


#endif
