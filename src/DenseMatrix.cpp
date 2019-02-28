/// \file DenseMatrix.cpp
/// Implementation of a DENSE matrix as
/// an Vector of DenseVector.

#include <complex>
#include <algorithm>

#include <DenseVector.h>
#include <DenseMatrix.h>
#include <Exceptions.h>
#include <Functors.h>
#include <Utility.h>

namespace CppNoddy {

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix() : m_nr(0), m_nc(0)
  {}

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix(const std::size_t& rows,
                                  const std::size_t& cols,
                                  const _Type& fill) :
    Sequential_Matrix_base<_Type>(),
    m_nr(rows),
    m_nc(cols) {
    // make a row
    const DenseVector<_Type> row(cols, fill);
    // reserve the space
    m_matrix.reserve(rows);
    for(std::size_t i = 0; i < rows; ++i) {
      // push require number of rows into the 'matrix'
      m_matrix.push_back(row);
    }
  }

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix(const std::size_t& rows,
                                  const std::size_t& cols,
                                  const _Type* p) :
    Sequential_Matrix_base<_Type>(),
    m_nr(rows),
    m_nc(cols) {
    m_matrix.reserve(rows);
    for(std::size_t i = 0; i < rows; ++i) {
      m_matrix.push_back(DenseVector<_Type>(cols, &p[ i * cols ]));
    }
  }

  template <typename _Type>
  DenseMatrix<_Type>::DenseMatrix(const DenseMatrix<_Type>& source) :
    Sequential_Matrix_base<_Type>() {
    *this = source;
  }

  template <typename _Type>
  DenseMatrix<_Type>::~DenseMatrix()
  {}

  template <typename _Type>
  inline DenseMatrix<_Type>& DenseMatrix<_Type>::operator=(const DenseMatrix<_Type>& source) {
    if(this == &source)
      return * this;
    m_matrix = source.m_matrix;
    m_nr = source.m_nr;
    m_nc = source.m_nc;
    return *this;
  }

  template <typename _Type>
  std::size_t DenseMatrix<_Type>::nelts() const {
    return m_nr * m_nc;
  }

  template <typename _Type>
  void DenseMatrix<_Type>::add(const DenseMatrix<_Type>& B) {
#ifdef PARANOID
    // check number of columns at least match
    if((B.nrows() != m_nr) || (B.ncols() != m_nc)) {
      std::string problem("The DenseMatrix.add has a geometry error.\n");
      throw ExceptionGeom(problem, m_nr, m_nc, B.nrows(), B.ncols());
    }
#endif
    std::transform(m_matrix.begin(), m_matrix.end(), B.m_matrix.begin(), m_matrix.begin(), std::plus< DenseVector<_Type> >());
  }

  template <typename _Type>
  void DenseMatrix<_Type>::sub(const DenseMatrix<_Type>& B) {
#ifdef PARANOID
    // check number of columns at least match
    if((B.nrows() != m_nr) || (B.ncols() != m_nc)) {
      std::string problem("The DenseMatrix.sub has a geometry error.\n");
      throw ExceptionGeom(problem, m_nr, m_nc, B.nrows(), B.ncols());
    }
#endif
    std::transform(m_matrix.begin(), m_matrix.end(), B.m_matrix.begin(), m_matrix.begin(), std::minus< DenseVector<_Type> >());
  }

  template <typename _Type>
  void DenseMatrix<_Type>::scale(const _Type& mult) {
    std::transform(m_matrix.begin(), m_matrix.end(), m_matrix.begin(), scale_functor< DenseVector<_Type>, _Type >(mult));
  }

  template <typename _Type>
  void DenseMatrix<_Type>::transpose() {
    if(nrows() == ncols()) {
      // square matrix needs no temp object
      // loop through upper half diagonal of the matrix
      for(std::size_t i = 0; i < nrows(); ++i) {
        for(std::size_t j = i + 1; j < ncols(); ++j) {
          // swap elements
          std::swap(m_matrix[ i ][ j ], m_matrix[ j ][ i ]);
        }
      }
    } else {
      std::vector< DenseVector<_Type> > temp;
      temp.resize(m_nc);
      for(std::size_t row = 0; row < m_nc; ++row) {
        temp[ row ].resize(m_nr);
      }
      for(std::size_t i = 0; i < nrows(); ++i) {
        for(std::size_t j = 0; j < ncols(); ++j) {
          temp[ j ][ i ] = m_matrix[ i ][ j ];
        }
      }
      m_matrix = temp;
      std::swap(m_nr, m_nc);
    }
  }

  template <typename _Type>
  DenseVector<_Type> DenseMatrix<_Type>::multiply(const DenseVector<_Type>& X) const {
#ifdef PARANOID
    // check number of columns at least match
    if(X.size() != m_nc) {
      std::string problem("The DenseMatrix.multiply has a geometry error.\n");
      throw ExceptionGeom(problem, m_nr, m_nc, X.size(), 1);
    }
#endif
    DenseVector<_Type> temp;
    temp.reserve(m_nr);
    for(std::size_t row = 0; row < m_nr; ++row) {
      temp.push_back(Utility::dot(m_matrix[ row ], X));
    }
    return temp;
  }

  template <typename _Type>
  DenseMatrix<_Type> DenseMatrix<_Type>::multiply(const DenseMatrix<_Type>& B) const {
#ifdef PARANOID
    // check number of columns at least match
    if(B.nrows() != m_nc) {
      std::string problem("The DenseMatrix.multiply has a geometry error.\n");
      throw ExceptionGeom(problem, m_nr, m_nc, B.nrows(), B.ncols());
    }
#endif
    // temporary object for the result
    DenseMatrix<_Type> C(nrows(), B.ncols(), 0.0);
    // loops thru the columns in the B matrix
    for(std::size_t col_in_B = 0; col_in_B < B.ncols(); ++col_in_B) {
      // set the column in the result to be the matrix-vector
      // product of (*this).multiply( column in B )
      C.set_col(col_in_B, multiply(B.get_col(col_in_B)));
    }
    return C;
  }


  template <typename _Type>
  typename std::vector<DenseVector<_Type> >::iterator
  DenseMatrix<_Type>::max_in_col(const std::size_t& col,
                                 row_iter row_min, row_iter row_max) const {
    row_iter index(row_min);
    double maxelt(std::abs(*(row_min -> begin())));
    for(row_iter row = row_min + 1; row != row_max ; ++row) {
      const double elt(std::abs(*(row -> begin() + col)));
      if(elt >= maxelt) {
        maxelt = elt;
        index = row;
      }
    }
    return index;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::one_norm() const {
    double max(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      max = std::max(max, m_matrix[ row ].one_norm());
    }
    return max;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::two_norm() const {
    double max(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      max = std::max(max, m_matrix[ row ].two_norm());
    }
    return max;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::inf_norm() const {
    double max(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      max = std::max(max, m_matrix[ row ].inf_norm());
    }
    return max;
  }

  template <typename _Type>
  double DenseMatrix<_Type>::frob_norm() const {
    double sum(0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      sum += m_matrix[ row ].two_norm();
    }
    return sum;
  }


  template <typename _Type>
  void DenseMatrix<_Type>::dump() const {
    std::cout << "DENSE mtx size = " << nrows() << " x " << ncols() << "; \n";
    std::cout.precision(3);
    std::cout << std::fixed;
    std::cout.setf(std::ios::showpoint);
    std::cout.setf(std::ios::showpos);
    //std::cout.setf( std::ios::scientific );
    std::cout << "- start matrix \n";
    for(std::size_t i = 0; i < nrows(); ++i) {
      std::cout << " row " << i << " =  ";
      for(std::size_t j = 0; j < ncols(); ++j) {
        std::cout << m_matrix[ i ][ j ] << ", ";
      }
      std::cout << "\n";
    }
    std::cout << "- end matrix \n";
  }

  template <typename _Type>
  void DenseMatrix<_Type>::set_col(const std::size_t& col, const DenseVector<_Type>& X) {
    for(std::size_t row = 0; row < m_nr; ++row) {
      m_matrix[ row ][ col ] = X[ row ];
    }
  }

  template <typename _Type>
  DenseVector<_Type> DenseMatrix<_Type>::get_col(const std::size_t& col) const {
    DenseVector<_Type> X(m_nr, 0.0);
    for(std::size_t row = 0; row < m_nr; ++row) {
      X[ row ] = m_matrix[ row ][ col ];
    }
    return X;
  }

  template <>
  void DenseMatrix<std::complex<double> >::matrix_to_vector(DenseVector<double> &p, const std::size_t &padding) const {
    p.reserve(2 * m_nr * m_nc);
    for(std::size_t row = 0; row < m_nr; ++row) {
      for(std::size_t col = 0; col < m_nc; ++col) {
        p.push_back(m_matrix[ row ][ col ].real());
        p.push_back(m_matrix[ row ][ col ].imag());
      }
      for(std::size_t col = 0; col < padding; ++col) {
        p.push_back(0.0);
        p.push_back(0.0);
      }
    }
  }


  template <>
  void DenseMatrix<double>::matrix_to_vector(DenseVector<double> &p, const std::size_t &padding) const {
    p.reserve(m_nr * m_nc);
    for(std::size_t row = 0; row < m_nr; ++row) {
      for(std::size_t col = 0; col < m_nc; ++col) {
        p.push_back(m_matrix[ row ][ col ]);
      }
      for(std::size_t col = 0; col < padding; ++col) {
        p.push_back(0.0);
      }
    }
  }

  template <>
  DenseVector<double> DenseMatrix<double>::matrix_to_vector(const std::size_t &padding) const {
    DenseVector<double> V;
    V.reserve(m_nr * m_nc);
    for(std::size_t row = 0; row < m_nr; ++row) {
      for(std::size_t col = 0; col < m_nc; ++col) {
        V.push_back(m_matrix[ row ][ col ]);
      }
      for(std::size_t col = 0; col < padding; ++col) {
        V.push_back(0.0);
      }
    }
    return V;
  }

  template <>
  DenseVector<double> DenseMatrix<std::complex<double> >::matrix_to_vector(const std::size_t &padding) const {
    DenseVector<double> V;
    V.reserve(2 * m_nr * m_nc);
    for(std::size_t row = 0; row < m_nr; ++row) {
      for(std::size_t col = 0; col < m_nc; ++col) {
        V.push_back(m_matrix[ row ][ col ].real());
        V.push_back(m_matrix[ row ][ col ].imag());
      }
      for(std::size_t col = 0; col < padding; ++col) {
        V.push_back(0.0);
        V.push_back(0.0);
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
