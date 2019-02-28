/// \file SparseMatrix.cpp
/// Implementation of a SPARSE matrix as
/// an STL Vector of SparseVectors, inheriting from Matrix_base.

#include <string>
#include <complex>

#include <SparseMatrix.h>
#include <Exceptions.h>


namespace CppNoddy {

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix(const std::size_t& rows, const std::size_t& cols)
    : m_nr(rows), m_nc(cols) {
    m_matrix.reserve(m_nr);
    SparseVector<_Type> sparse_row(m_nc);
    for(std::size_t i = 0; i < m_nr; ++i) {
      m_matrix.push_back(sparse_row);
    }
  }

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix(const SparseMatrix<_Type>& source, const std::vector<std::size_t>& source_rows) :
    m_nr(source.nrows()), m_nc(source.ncols()) {
    m_matrix.reserve(m_nr);
    for(std::size_t i = 0; i < m_nr; ++i) {
      m_matrix.push_back(source.get_row(source_rows[i]));
    }
  }

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix(const SparseMatrix<_Type>& source) {
    *this = source;
  }

  template <typename _Type>
  inline SparseMatrix<_Type>& SparseMatrix<_Type>::operator=(const SparseMatrix<_Type>& source) {
    if(this == &source)
      return * this;
    m_matrix = source.m_matrix;
    m_nr = source.m_nr;
    m_nc = source.m_nc;
    return *this;
  }

  template <typename _Type>
  std::size_t SparseMatrix<_Type>::nelts() const {
    std::size_t num_of_elts(0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      num_of_elts += m_matrix[ row ].nelts();
    }
    return num_of_elts;
  }

  template <typename _Type>
  std::size_t SparseMatrix<_Type>::max_in_col(const std::size_t& col,
      const std::size_t& row_min, const std::size_t& row_max) const {
    double maxelt(0.0);
    // return outside of the array as default
    std::size_t index = nrows();
    for(std::size_t row = row_min; row < row_max ; ++row) {
      // only bother looking up entries of rows with a first
      // element in a column less than the one we're looking at
      if(m_matrix[ row ].begin() -> first <= col) {
        const double elt(std::abs(m_matrix[ row ].get(col)));
        if(elt >= maxelt) {
          maxelt = elt;
          index = row;
        }
      }
    }
    return index;
  }

  template <typename _Type>
  void SparseMatrix<_Type>::scale(const _Type& mult) {
    for(std::size_t row = 0; row < m_nr; ++row) {
      m_matrix[ row ] *= mult;
    }
  }
  
  template <typename _Type>
  double SparseMatrix<_Type>::one_norm() const {
    double max(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      max = std::max(max, m_matrix[ row ].one_norm());
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::two_norm() const {
    double max(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      max = std::max(max, m_matrix[ row ].two_norm());
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::inf_norm() const {
    double max(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      max = std::max(max, m_matrix[ row ].inf_norm());
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::frob_norm() const {
    double sum(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      sum += m_matrix[ row ].two_norm();
    }
    return sum;
  }

  
#if defined(PETSC_Z)
  template <>
  void SparseMatrix<std::complex<double> >::get_row_petsc(PetscInt row,
                  PetscScalar* storage, PetscInt* cols) {
    citer pos;
    std::size_t i(0);
    //
    // matrix could be singular with an empty row for the mass matrix
    // of a generalised eigenvalue problem
    if(m_matrix[row].nelts() > 0) {
      // start at the begining of this row
      pos = m_matrix[ row ].begin();
      do {
        // for each non-zero elt in the row
        PetscScalar elt;
        elt = std::real(pos -> second) + PETSC_i * std::imag(pos -> second);
        int col(pos -> first);
        storage[ i ] = elt;
        // +1 to return FORTRAN indexing
        cols[ i ] = col;
        ++pos;
        ++i;
      } while(pos != m_matrix[ row ].end());
    }
  }
#endif

#if defined(PETSC_D)
  template <>
  void SparseMatrix<std::complex<double> >::get_row_petsc(PetscInt row, PetscScalar* storage, PetscInt* cols) {
    std::string problem;
    problem = "The SparseMatrix::get_row_petsc method was called for a SparseMatrix<D_complex>\n";
    problem += "even though PETSC_ARCH is currently pointing to a double version of the library.\n";
    throw ExceptionExternal(problem);
  }
#endif


#if defined(PETSC_D)
  template <>
  void SparseMatrix<double >::get_row_petsc(PetscInt row, PetscScalar* storage, PetscInt* cols) {
    citer pos;
    std::size_t i(0);
    //
    // matrix could be singular with an empty row for the mass matrix
    // of a generalised eigenvalue problem
    if(m_matrix[row].nelts() > 0) {
      // start at the begining of this row
      pos = m_matrix[ row ].begin();
      do {
        // for each non-zero elt in the row
        PetscScalar elt;
        elt = pos -> second;
        int col(pos -> first);
        storage[ i ] = elt;
        // +1 to return FORTRAN indexing
        cols[ i ] = col;
        ++pos;
        ++i;
      } while(pos != m_matrix[ row ].end());
    }
  }
#endif

#if defined(PETSC_Z)
  template <>
  void SparseMatrix<double >::get_row_petsc(PetscInt row, PetscScalar* storage, PetscInt* cols) {
    std::string problem;
    problem = "The SparseMatrix::get_row_petsc method was called for a SparseMatrix<double>\n";
    problem += "even though PETSC_ARCH is currently pointing to a complex version of the library.\n";
    throw ExceptionExternal(problem);
  }
#endif


  template <typename _Type>
  void SparseMatrix<_Type>::dump() const {
    std::cout << "SPARSE mtx size = " << m_nr << " x  sparse \n";
    std::cout.precision(4);
    std::cout << "- start matrix \n";
    for(std::size_t row = 0; row < m_nr; ++row) {
      std::cout << " row " << row << " :  ";
      m_matrix[ row ].dump();
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
