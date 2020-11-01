#pragma once

#include <SparseVector.h>
#include <DenseVector.h>
#include <Sequential_Matrix_base.h>
#include <BandedMatrix.h>
#include <DenseMatrix.h>
#include <ctime>

namespace Utils_Fill {


  /// initialise RNG
  void time_seed() {
    srand((unsigned) std::time(0));
  }

  
  /// Fill diagonal with unit values
  /// \param A The dense matrix to be filled
  template <typename _Type>
  void fill_identity(CppNoddy::Sequential_Matrix_base<_Type>& A ) {
    for (std::size_t row = 0; row < A.nrows(); ++row) {
      A(row, row) = 1.0;
    }
  }
  
  /// Fill a diagonal band of a matrix
  /// \param A The matrix to be used
  /// \param offset The offset of the band from the main diagonal
  ///  e.g. 0 = main diagional, -1 = first sub-diagonal
  /// \param value The value to be written to the band elements
  template <typename _Type>
  void fill_band(CppNoddy::Sequential_Matrix_base<_Type>& A, const int& offset, const _Type& value) {
    for (std::size_t row = 0; row < A.nrows(); ++row) {
      if ((row + offset < A.ncols()) && (row + offset >= 0)) {
        A(row, row + offset) = value;
      }
    }
  }
  
  /// Set all elements of a DENSE vector
  /// \param X The DENSE vector to be filled
  /// \param value The value to be placed in each element of the vector
  template <typename _Type>
  void fill(CppNoddy::DenseVector<_Type>& X, const _Type& value) {
    for(std::size_t i = 0; i < X.size(); ++i) {
      X[ i ] = value;
    }
  }
  
  
  void fill_random(CppNoddy::SparseVector<double>& V, const unsigned& num_of_elts) {
    do {
      double index = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      index *= V.size();
      double x = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      V[(unsigned) index ] = x;
    } while(V.nelts() < num_of_elts);
  }
  
  void fill_random(CppNoddy::SparseVector<std::complex<double> >& V, const unsigned& num_of_elts) {
    do {
      double index = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      index *= V.size();
      double x = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      double y = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      V[(unsigned) index ] = std::complex<double>(x, y);
    } while(V.nelts() < num_of_elts);
  }
  
  void fill_random(CppNoddy::DenseVector<double>& V) {
    for(unsigned i = 0; i < V.size(); ++i) {
      double x = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      V[ i ] = x;
    }
  }

  void fill_random(CppNoddy::DenseVector<std::complex<double> >& V) {
    for(unsigned i = 0; i < V.size(); ++i) {
      double x = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      double y = (double) rand() /
        ((double) RAND_MAX + (double) 1) ;
      V[ i ] = std::complex<double>(x, y);
    }
  }
  
  void fill_random(CppNoddy::DenseMatrix<double>& A) {
    CppNoddy::DenseVector<double> temp(A.ncols(), 0.0);
    for(std::size_t row = 0; row < A.nrows(); ++row) {
      fill_random(temp);
      A[ row ] = temp;
    }
  }
  
  void fill_random(CppNoddy::BandedMatrix<double>& A) {
    for(std::size_t row = 0; row < A.nrows(); ++row) {
      for(std::size_t col = std::max((int)(row - A.noffdiag()), 0);
          (int) col <= std::min((int)(row + A.noffdiag()), (int) A.ncols()); ++col) {
        double x = (double) rand() / ((double) RAND_MAX + (double) 1) ;
        A(row, col) = x;
      }
    }
  }
  
}
