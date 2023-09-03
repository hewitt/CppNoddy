/// \file BandedMatrix.cpp
/// Implementation of the matrix class that constructs a BANDED matrix as
/// a DenseVector.

#include <string>
#include <memory>
#include <algorithm>

#include <DenseVector.h>
#include <BandedMatrix.h>
#include <Exceptions.h>

namespace CppNoddy {

  template <typename _Type>
  BandedMatrix<_Type>::BandedMatrix(const std::size_t& rows,
                                    const std::size_t& upper_offdiag_bands,
                                    const _Type& fill) :
    Sequential_Matrix_base<_Type>(),
    m_N(rows),
    m_L(upper_offdiag_bands) {
    // we'll assume that the upper & lower offdiags are of equal size.
    // logical m_storage is the main diagonal (=1) + upper offdiagonal (=m_L) + lower offdiagonal (=m_L)
    // = 2 * m_L + 1.
    // However, allowing for row interchanging during pivotting leads to additional
    // padding of m_L, to make a total of 3*m_L+1.
    m_storage = DenseVector<_Type>(m_N * (3 * m_L + 1), 0.0);
  }

  template <typename _Type>
  BandedMatrix<_Type>::BandedMatrix(const BandedMatrix<_Type>& source) :
    Sequential_Matrix_base<_Type>() {
    *this = source;
  }

  template <typename _Type>
  inline BandedMatrix<_Type>& BandedMatrix<_Type>::operator=(const BandedMatrix<_Type>& source) {
    if(this == &source)
      return * this;
    m_storage = source.m_storage;
    m_N = source.m_N;
    m_L = source.m_L;
    return *this;
  }
    
  template <typename _Type>
  std::size_t BandedMatrix<_Type>::nelts() const {
    return m_storage.size();
  }

  template <typename _Type>
  void BandedMatrix<_Type>::scale(const _Type& mult) {
    m_storage *= mult;
  }

  template <typename _Type>
  DenseVector<_Type> BandedMatrix<_Type>::multiply(const DenseVector<_Type>& X) const {
    std::string problem;
    problem = " The multiply method of the BandedMatrix class has \n";
    problem += " not been implemented. \n";
    throw ExceptionRuntime(problem);
    return m_storage; // dummy return
  }


  template <typename _Type>
  void BandedMatrix<_Type>::transpose() {
    for(std::size_t row = 0; row < m_N; ++row) {
      for(std::size_t col = std::max(0, int(row) - int(m_L)); col < row; ++col) {
        // swap elements
        std::swap(operator()(row, col), operator()(col, row));
      }
    }
  }


  template <typename _Type>
  void BandedMatrix<_Type>::row_swap(const std::size_t& row1, const std::size_t& row2) {
    // WE MUST HAVE row2 > row1 & row1 must have no elements left of the diagonal!
    // only of any use to native Gaussian elimination solver
    /// \todo MAKE ME PRIVATE OR ELSE!
    for(std::size_t col = row1; col <= std::min(row1 + 2 * m_L, m_N - 1); ++col) {
      std::swap(operator()(row1, col), operator()(row2, col));
    }
  }


  template <typename _Type>
  double BandedMatrix<_Type>::one_norm() const {
    double max(0.0);
    for(unsigned row = 0; row < m_N; ++row) {
      // copy the row into a temp vector
      DenseVector<_Type> temp(3 * m_L + 1, &m_storage[ row * (3 * m_L + 1) ]);
      max = std::max(max, temp.one_norm());
    }
    return max;
  }


  template <typename _Type>
  double BandedMatrix<_Type>::two_norm() const {
    double max(0.0);
    for(unsigned row = 0; row < m_N; ++row) {
      // copy the row into a temp vector
      DenseVector<_Type> temp(3 * m_L + 1, &m_storage[ row * (3 * m_L + 1) ]);
      max = std::max(max, temp.two_norm());
    }
    return max;
  }


  template <typename _Type>
  double BandedMatrix<_Type>::inf_norm() const {
    double max(0.0);
    for(unsigned row = 0; row < m_N; ++row) {
      // copy the row into a temp vector
      DenseVector<_Type> temp(3 * m_L + 1, &m_storage[ row * (3 * m_L + 1) ]);
      max = std::max(max, temp.inf_norm());
    }
    return max;
  }


  template <typename _Type>
  double BandedMatrix<_Type>::frob_norm() const {
    double sum(0.0);
    for(unsigned row = 0; row < m_N; ++row) {
      // copy the row into a temp vector
      DenseVector<_Type> temp(3 * m_L + 1, &m_storage[ row * (3 * m_L + 1) ]);
      sum += temp.two_norm();
    }
    return sum;
  }

  template <typename _Type>
  void BandedMatrix<_Type>::dump() const {
    std::cout << "BANDED mtx size = " << nrows() << " x " << ncols()
              << " with INPUT no. of bands = " << 2 * m_L + 1
              << " and total bands STORED (for pivotting) = " << 3 * m_L + 1 << "\n";
    std::cout.precision(4);
    std::cout << "- start matrix \n";
    for(std::size_t i = 0; i < nrows(); ++i) {
      std::cout << " row " << i << " =  ";
      for(std::size_t j = 0; j < nrows(); ++j) {
        if(i == j)
          std::cout << "*";
        if(((int) j < (int) i - (int) m_L) || ((int) j > (int) i + (int)(2*m_L))) {
          std::cout << 0 << ", ";
        } else {
          std::cout << operator()(i, j) << ", ";
        }
      }
      std::cout << "\n";
    }
    std::cout << "- end matrix \n";
  }

  template <>
  double* BandedMatrix<double>::base() {
    return &(m_storage[0]);
  }

  template <>
  double* BandedMatrix<std::complex<double> >::base() {
    //return &( m_storage[0].real() );
    return &reinterpret_cast<double(&)[2]>(m_storage[0])[0];
  }

  // the templated versions we require are:
  template class BandedMatrix<double>
  ;
  template class BandedMatrix<std::complex<double> >
  ;

} // end namespace
