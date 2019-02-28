
/// \file BandedLinearSystem.cpp
/// Implementation for the LinearSystem class

#include <vector>
#include <set>
#include <cassert>

#include <FortranLAPACK.h>
#include <BandedLinearSystem.h>
#include <Exceptions.h>
#include <Utility.h>
#include <Types.h>

namespace CppNoddy {

  template <typename _Type>
  BandedLinearSystem<_Type>::BandedLinearSystem
  (BandedMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which) :
    m_detSign(0), m_monitorDet(false) {
    m_pA = Aptr;
    m_pB = Bptr;
    m_version = which;
    if((m_version != "lapack") && (m_version != "native")) {
      std::string problem;
      problem = "The BandedLinearSystem has been instantiated with an unrecognised\n";
      problem += "request for a solver type. Options are 'native' or 'lapack'. \n";
      throw ExceptionRuntime(problem);
    }
  }
  
  
  template <typename _Type>
  void BandedLinearSystem<_Type>::solve() {
    if("lapack" == m_version) {
      solve_lapack();
    } else { // we catch incorrect m_version choices in the ctor
      solve_native();
    }
  }

  // only the specialised solvers below are available
  template <typename _Type>
  void BandedLinearSystem<_Type>::solve_lapack() {
    std::string problem;
    problem = "The solve method for a BandedLinearSystem has not been implemented\n";
    problem += "for the element type used here. \n";
    throw ExceptionExternal(problem);
  }

  // lapack LU solver for double-element banded matrices
  template <>
  void BandedLinearSystem<double>::solve_lapack() {
#ifndef LAPACK
    std::string problem = "The BandedLinearSystem<double>::solve_lapack method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built and so LAPACK support is not available.";
    throw ExceptionExternal(problem);
#else
    const std::size_t N = m_pA -> nrows();
    const std::size_t K = 1;
    const std::size_t L = m_pA -> noffdiag();
    const std::size_t LDAB = 3 * L + 1;
    // pivot storage
    std::vector<int> pivots(N);
    // warning info
    int info(0);
    LAPACK_DGBSV(N, L, L, K, m_pA -> base(), LDAB, &pivots[ 0 ], &(m_pB -> operator[](0)), N, info);
    m_pivots = pivots;
    if(0 != info) {
      std::string problem(" The BandedLinearSystem::solve_lapack method has detected a failure. \n");
      throw ExceptionExternal(problem, info);
    }
    // compute the determinant if asked
    if(m_monitorDet) {
      m_detSign = signature(pivots);
      // product of the diagonals
      for(std::size_t i = 0; i < N; ++i) {
        if((*m_pA)(i, i) < 0.0) {
          m_detSign *= -1;
        }
      }
    }
#endif
  }

  // lapack LU re-solver for double-element banded matrices
  template <>
  void BandedLinearSystem<double>::re_solve_lapack() {
#ifndef LAPACK
    std::string problem = "The BandedLinearSystem<double>::re_solve_lapack method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built and so LAPACK support is not available.";
    throw ExceptionExternal(problem);
#else
    const std::size_t N = m_pA -> nrows();
    const std::size_t K = 1;
    const std::size_t L = m_pA -> noffdiag();
    const std::size_t LDAB = 3 * L + 1;
    // warning info
    int info(0);
    // LAPACK_DGBSV( N, L, L, K, m_pA -> base(), LDAB, &pivots[ 0 ], &( m_pB -> operator[] ( 0 ) ), N, info );
    LAPACK_DGBTRS((char*) "N", N, L, L, K, m_pA -> base(), LDAB, &m_pivots[ 0 ], &(m_pB -> operator[](0)), N, info);
    if(0 != info) {
      std::string problem(" The BandedLinearSystem::re_solve_lapack method has detected a failure. \n");
      throw ExceptionExternal(problem, info);
    }
#endif
  }

  
  // lapack LU solver for complex-element banded matrices
  template <>
  void BandedLinearSystem<D_complex>::solve_lapack() {
#ifndef LAPACK
    std::string problem = "The BandedLinearSystem<D_complex>::solve_lapack method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built and so LAPACK support is not available.";
    throw ExceptionExternal(problem);
#else
    const std::size_t N = m_pA -> nrows();
    const std::size_t K = 1;
    const std::size_t L = m_pA -> noffdiag();
    const std::size_t LDAB = 3 * L + 1;
    // pivot storage
    std::vector<int> pivots(N);
    // warning info
    int info(0);
    // LAPACK_ZGBSV( N, L, L, K, m_pA -> base(), LDAB, &pivots[ 0 ], &( m_pB -> operator[]( 0 ).real() ), N, info );
    LAPACK_ZGBSV(N, L, L, K, m_pA -> base(), LDAB, &pivots[ 0 ], &reinterpret_cast<double(&)[2]>(m_pB -> operator[](0))[0], N, info);
    if(0 != info) {
      std::string problem(" The BandedLinearSystem::solve_lapack method has detected a failure. \n");
      throw ExceptionExternal(problem, info);
    }
#endif
  }


  template <typename _Type>
  void BandedLinearSystem<_Type>::solve_native() {
    // determinant sign correction from pivotting
    int sign(1);
    // number of offdiagonal elts STORED above the main diagonal
    const unsigned L = m_pA -> noffdiag();
    const std::size_t N = m_pA -> nrows();
    // step through rows
    for(std::size_t l = 0 ; l < N - 1 ; ++l) {
      _Type diag_entry = (*m_pA)(l, l);
      // eliminate all entries below
      // this is the maximum number of non-zero elts below the diagonal
      // +1 because of usual zero indexing
      const std::size_t row_end = std::min(N, l + 1 * L + 1);
      const std::size_t col_end = std::min(N, l + 2 * L + 1);

      // assume the diagonal elt is currently the maximum one
      std::size_t max_row = l;
      // partial pivot by looking at this COLUMN for bigger values
      for(std::size_t row = l + 1; row < row_end; ++row) {
        if(std::abs((*m_pA)(row, l)) > std::abs(diag_entry)) {
          diag_entry = (*m_pA)(row, l);
          max_row = row;
        }
      }
      // swap rows if suitable maximum is found
      if(max_row != l) {
        sign *= -1;
        m_pA -> row_swap(l, max_row);
        // swap RHS
        m_pB -> swap(l, max_row);
      }
      // eliminate
      for(std::size_t row = l + 1 ; row < row_end; ++row) {
        // work out mulitplier
        const _Type mult = (*m_pA)(row, l) / diag_entry;
        // how to subtract two rows
        for(std::size_t col = l; col < col_end; ++col) {
          (*m_pA)(row, col) -= (*m_pA)(l, col) * mult;
        }
        // do the same subtraction of the RHS
        m_pB -> operator[](row) -= mult * m_pB -> operator[](l);
      } // close row-loop
    }  // close l-loop

    if(m_monitorDet) {
      m_detSign = sign;
      for(std::size_t i = 0; i < N; ++i) {
        if(lt(m_pA -> operator()(i, i))) {
          m_detSign *= -1;
        }
      }
    }

    // backsubstitute upper triangular matrix
    backsub();
  }

  template <typename _Type>
  bool BandedLinearSystem<_Type>::lt(_Type value) const {
    std::string problem;
    problem = "You've turned on monitoring of the sign of the determinant for a \n";
    problem += "BandedLinearSystem whose elements are not of type <double>.\n";
    throw ExceptionRuntime(problem);
  }

  template <>
  bool BandedLinearSystem<double>::lt(double value) const {
    return (value < 0);
  }

  template <typename _Type>
  int BandedLinearSystem<_Type>::signature(const std::vector<int> &pivots) const {
    int sign(1);
    for(std::size_t i = 0; i < pivots.size(); ++i) {
      if(pivots[ i ] - 1 != int(i)) {
        sign *= -1;
      }
    }
    return sign;
  }

  template <typename _Type>
  void BandedLinearSystem<_Type >::backsub() const {
    const std::size_t L = m_pA -> noffdiag();
    const std::size_t N = m_pB -> size();
    (*m_pB)[ N - 1 ] /= (*m_pA)(N - 1, N - 1);
    // Note the unusual row termination condition.
    // We can't do just "row >= 0" with std::size_t
    for(std::size_t row = N - 2; (int) row >= 0; --row) {
      _Type sum(0.0);
      for(std::size_t col = row + 1; col < std::min(N, row + 2 * L + 1); ++col) {
        sum += (*m_pA)(row, col) * (*m_pB)[ col ];
      }
      (*m_pB)[ row ] -= sum;
      (*m_pB)[ row ] /= (*m_pA)(row, row);
    }
  }

  template<typename _Type>
  int BandedLinearSystem<_Type>::get_det_sign() const {
    return m_detSign;
  }

  template<typename _Type>
  void BandedLinearSystem<_Type>::set_monitor_det(bool flag) {
    m_monitorDet = flag;
  }

  

  template class BandedLinearSystem<D_complex>
  ;
  template class BandedLinearSystem<double>
  ;

} // end namespace
