/// \file DenseLinearEigenSystem.cpp
/// Implementation for the DenseLinearEigenSystem class
/// This class links to LAPACK to perform the solver phase.

#include <vector>
#include <set>

#include <FortranLAPACK.h>
#include <FortranData.h>
#include <DenseLinearEigenSystem.h>
#include <Exceptions.h>
#include <Types.h>

namespace CppNoddy {

  template <typename _Type>
  DenseLinearEigenSystem<_Type>::DenseLinearEigenSystem(DenseMatrix<_Type>* Aptr, DenseMatrix<_Type>* Bptr) :
    LinearEigenSystem_base() {
    m_pA = Aptr;
    m_pB = Bptr;
  }

  template <typename _Type>
  DenseLinearEigenSystem<_Type>::~DenseLinearEigenSystem()
  {}

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::eigensolve() {
    // reset any tagged eigenvalues
    m_tagged_indices.clear();
    //
    if(m_calc_eigenvectors) {
      eigensolve_lapack_with_vectors();
    } else {
      eigensolve_lapack_without_vectors();
    }
  }


  template <>
  void DenseLinearEigenSystem<double>::eigensolve_lapack_without_vectors() {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal(problem);
#else

    std::size_t N = m_pA -> nrows();
    // Cache issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding(0);
    if((N % 2 == 0) && (N > 127)) {
      padding = 1;
    }
#ifdef PARANOID
    if((m_pA -> nrows() != m_pB -> nrows()) ||
        (m_pA -> ncols() != m_pB -> ncols())) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure \n");
      throw ExceptionGeom(problem, m_pA -> nrows(), m_pA -> ncols(),
                          m_pB -> nrows(), m_pB -> ncols());
    }
#endif
    FortranData Af(*m_pA, true, padding);
    FortranData Bf(*m_pB, true, padding);
    // eigenvalue storage
    DenseVector<double> alpha_r(N, 0.0);
    DenseVector<double> alpha_i(N, 0.0);
    DenseVector<double> beta(N, 0.0);
    // eigenvector storage
    DenseVector<double> vec_left(1, 0.0);
    DenseVector<double> vec_right(1, 0.0);
    // some workspace for the LAPACK routine
    DenseVector<double> work(1, 0.0);
    int info(0);
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_DGGEV((char*) "N", (char*) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], -1, info);
    int required_workspace = (int) work[ 0 ];
#ifdef DEBUG

    std::cout << " [DEBUG] DenseLinearEigenSystem::eigensolve_lapack_without_vectors is requesting \n";
    std::cout << " [DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize(required_workspace);
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_DGGEV((char*) "N", (char*) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], required_workspace, info);
    if(0 != info) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure. \n");
      throw ExceptionExternal(problem, info);
    }
    // create a complex eigenvalue vector
    m_eigenvalues_alpha = DenseVector<D_complex>(N, 0.0);
    for(std::size_t i = 0; i < N; ++i) {
      const D_complex eye(0.0, 1.0);
      m_eigenvalues_alpha[ i ] = alpha_r[ i ] + alpha_i[ i ] * eye;
    }
    // set the eigenvalue member data
    m_eigenvalues_beta = beta;
#endif

  }

  template <>
  void DenseLinearEigenSystem< std::complex<double> >::eigensolve_lapack_without_vectors() {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal(problem);
#else

    std::size_t N = m_pA -> nrows();
    // Cache issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding(0);
    if((N % 2 == 0) && (N > 127)) {
      padding = 1;
    }
#ifdef PARANOID
    if((m_pA -> nrows() != m_pB -> nrows()) ||
        (m_pA -> ncols() != m_pB -> ncols())) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure. \n");
      throw ExceptionGeom(problem, m_pA -> nrows(), m_pA -> ncols(),
                          m_pB -> nrows(), m_pB -> ncols());
    }
#endif
    // transpose the input matrices so they are in column_major format
    FortranData Af(*m_pA, true, padding);
    FortranData Bf(*m_pB, true, padding);
    // eigenvalue storage
    DenseVector<double> alpha(2 * N, 0.0);
    DenseVector<double> beta(2 * N, 0.0);
    // eigenvector storage
    DenseVector<double> vec_left(2, 0.0);
    DenseVector<double> vec_right(2, 0.0);
    // new the eigenvector storage
    // some workspace for the LAPACK routine
    DenseVector<double> work(2, 0.0);
    DenseVector<double> rwork(8 * N, 0.0);
    int info(0);
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_ZGGEV((char*) "N", (char*) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], -1, &rwork[ 0 ], info);
    int required_workspace = int(work[ 0 ]);
#ifdef DEBUG

    std::cout << "[DEBUG] DenseLinearEigenSystem::eigensolve_lapack_without_vectors is requesting \n";
    std::cout << "[DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize(2 * required_workspace);
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_ZGGEV((char*) "N", (char*) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], required_workspace, &rwork[ 0 ], info);
    // error reporting
    if(0 != info) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure \n");
      throw ExceptionExternal(problem, info);
    }
    // create a complex eigenvalue vector for returning
    m_eigenvalues_alpha = DenseVector<D_complex>(N, 0.0);
    m_eigenvalues_beta = DenseVector<D_complex>(N, 0.0);
    {
      const D_complex eye(0.0, 1.0);
      for(std::size_t i = 0; i < N; ++i) {
        m_eigenvalues_alpha[ i ] = alpha[ 2 * i ] + alpha[ 2 * i + 1 ] * eye;
        m_eigenvalues_beta[ i ] = beta[ 2 * i ] + beta[ 2 * i + 1 ] * eye;
      }
    }
#endif

  }

  template <>
  void DenseLinearEigenSystem< double >::eigensolve_lapack_with_vectors() {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal(problem);
#else

    std::size_t N = m_pA -> nrows();
    // Cache contention issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding(0);
    if((N % 2 == 0) && (N > 127)) {
      padding = 1;
    }
#ifdef PARANOID
    if((m_pA -> nrows() != m_pB -> nrows()) ||
        (m_pA -> ncols() != m_pB -> ncols())) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure. \n");
      throw ExceptionGeom(problem, m_pA -> nrows(), m_pA -> ncols(),
                          m_pB -> nrows(), m_pB -> ncols());
    }
#endif
    // Convert to fortran data  incl. a transpose of the
    // input matrices so they are in column_major format then include padding
    FortranData Af(*m_pA, true, padding);
    FortranData Bf(*m_pB, true, padding);
    // eigenvalue storage
    DenseVector<double> alpha_r(N, 0.0);
    DenseVector<double> alpha_i(N, 0.0);
    DenseVector<double> beta(N, 0.0);
    // new the right eigenvector storage
    DenseVector<double> vec_left(1, 0.0);
    DenseVector<double> vec_right(N * N, 0.0);
    // some workspace for the LAPACK routine
    DenseVector<double> work(1, 0.0);
    // return integer for LAPACK
    int info(0);
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_DGGEV((char*) "N", (char*) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], -1, info);
    int required_workspace = 4 * int(work[ 0 ]);
#ifdef DEBUG

    std::cout << "[DEBUG] DenseLinearEigenSystem::eigensolve_lapack_with_vectors is requesting \n";
    std::cout << "[DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize(required_workspace);
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_DGGEV((char*) "N", (char*) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], required_workspace, info);
    if(0 != info) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure.\n");
      throw ExceptionExternal(problem, info);
    }
    // create a complex eigenvalue vector
    m_eigenvalues_alpha = DenseVector<D_complex>(N, 0.0);
    // complex eigenvector matrix
    m_all_eigenvectors = DenseMatrix<D_complex>(N, N, 0.0);
    // step through the eigenvalues
    for(std::size_t i = 0; i < N; ++i) {
      const D_complex eye(0.0, 1.0);
      // make the complex vector of alpha
      m_eigenvalues_alpha[ i ] = alpha_r[ i ] + alpha_i[ i ] * eye;
      if(std::abs(alpha_i[ i ]) > 0.0) {
        // eigenvector is complex
        for(std::size_t k = 0; k < N; ++k) {
          m_all_eigenvectors[ i ][ k ] = vec_right[ i * N + k ] + eye * vec_right[(i + 1) * N + k ];
          // store the conjugate too for completeness
          m_all_eigenvectors[ i + 1 ][ k ] = vec_right[ i * N + k ] - eye * vec_right[(i + 1) * N + k ];
        }
        ++i;
      } else { // eigenvector is real
        for(std::size_t k = 0; k < N; ++k) {
          m_all_eigenvectors(i, k) = vec_right[ i * N + k ];
        }
      }
    }
    // set the eigenvalue member data
    m_eigenvalues_beta = beta;
#endif

  }

  template <>
  void DenseLinearEigenSystem< std::complex<double> >::eigensolve_lapack_with_vectors() {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal(problem);
#else

    std::size_t N = m_pA -> nrows();
    // Cache issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding(0);
    if((N % 2 == 0) && (N > 127)) {
      padding = 1;
    }
#ifdef PARANOID
    if((m_pA -> nrows() != m_pB -> nrows()) ||
        (m_pA -> ncols() != m_pB -> ncols())) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure. \n");
      throw ExceptionGeom(problem, m_pA -> nrows(), m_pA -> ncols(),
                          m_pB -> nrows(), m_pB -> ncols());
    }
#endif
    // transpose the input matrices so they are in column_major format
    FortranData Af(*m_pA, true, padding);
    FortranData Bf(*m_pB, true, padding);
    // eigenvalue storage
    DenseVector<double> alpha(2 * N, 0.0);
    DenseVector<double> beta(2 * N, 0.0);
    // eigenvector storage
    DenseVector<double> vec_left(2, 0.0);
    DenseVector<double> vec_right(2 * N * N, 0.0);
    // some workspace for the LAPACK routine
    DenseVector<double> work(2, 0.0);
    DenseVector<double> rwork(8 * N, 0.0);
    int info(0);
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_ZGGEV((char*) "N", (char*) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], -1, &rwork[ 0 ], info);
    int required_workspace = int(work[ 0 ]);
#ifdef DEBUG

    std::cout << "[DEBUG] DenseLinearEigenSystem::eigensolve_lapack_with_vectors is requesting \n";
    std::cout << "[DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize(2 * required_workspace);
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_ZGGEV((char*) "N", (char*) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], required_workspace, &rwork[ 0 ], info);
    if(0 != info) {
      std::string problem("The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure.\n");
      throw ExceptionExternal(problem, info);
    }
    // create a complex eigenvalue vector
    m_eigenvalues_alpha = DenseVector<D_complex>(N, 0.0);
    m_eigenvalues_beta = DenseVector<D_complex>(N, 0.0);
    // complex eigenvector matrix
    m_all_eigenvectors = DenseMatrix<D_complex>(N, N, 0.0);
    // step through the eigenvalues
    for(std::size_t i = 0; i < N; ++i) {
      const D_complex eye(0.0, 1.0);
      m_eigenvalues_alpha[ i ] = alpha[ 2 * i ] + alpha[ 2 * i + 1 ] * eye;
      m_eigenvalues_beta[ i ] = beta[ 2 * i ] + beta[ 2 * i + 1 ] * eye;
      for(std::size_t j = 0; j < N; ++j) {
        m_all_eigenvectors(i, j) = vec_right[ 2 * i * N + 2 * j ] + vec_right[ 2 * i * N + 2 * j + 1 ] * eye;
      }
    }
#endif

  }

  template <typename _Type>
  DenseVector<D_complex> DenseLinearEigenSystem<_Type>::get_tagged_eigenvalues() const {
    if(m_tagged_indices.size() == 0) {
      std::string problem;
      problem = "In DenseLinearEigenSystem.get_tagged_eigenvalues() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime(problem);
    }
    // storage for the eigenvalues
    DenseVector<D_complex> evals;
    // loop through the tagged set
    for(iter p = m_tagged_indices.begin(); p != m_tagged_indices.end(); ++p) {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      // work out the complex eigenvalue associated with this index
      // and add it to the vector
      evals.push_back(m_eigenvalues_alpha[ j ] / m_eigenvalues_beta[ j ]);
    }
    // return the complex vector of eigenvalues
    return evals;
  }

  template <typename _Type>
  DenseMatrix<D_complex> DenseLinearEigenSystem<_Type>::get_tagged_eigenvectors() const {
    if(m_tagged_indices.size() == 0) {
      std::string problem;
      problem = "In DenseLinearEigenSystem.get_tagged_eigenvectors() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime(problem);
    }
    // order of the problem
    std::size_t N = m_eigenvalues_alpha.size();
    // eigenvector storage : size() eigenvectors each of length N
    DenseMatrix<D_complex> evecs(m_tagged_indices.size(), N, 0.0);
    std::size_t row = 0;
    // loop through the tagged set
    for(iter p = m_tagged_indices.begin(); p != m_tagged_indices.end(); ++p) {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      // put the eigenvector in the matrix
      evecs[ row ] = m_all_eigenvectors[ j ];
      // next row/eigenvector
      ++row;
    }
    return evecs;
  }

  // EIGENVALUE/VECTOR TAGGING

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_disc(const int &val, const double& radius) {
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_eigenvalues_alpha.size(); ++i) {
      // if the eigenvalue is in the disc centred at m_shift then include it
      if(std::abs(m_eigenvalues_alpha[ i ] - m_shift * m_eigenvalues_beta[ i ]) <
          std::abs(radius * m_eigenvalues_beta[ i ])) {
        if(val > 0) {
          // add it to our set of tagged eigenvalues
          m_tagged_indices.insert(m_tagged_indices.end(), i);
        } else {
          // remove it from the set if it exists
          m_tagged_indices.erase(i);
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_right(const int &val) {
    double real_value = m_shift.real();
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_eigenvalues_alpha.size(); ++i) {
      // if the eigenvalue is to the right of the shift position
      if((m_eigenvalues_alpha[ i ] * std::conj(m_eigenvalues_beta[ i ])).real() >
          std::pow(std::abs(m_eigenvalues_beta[ i ]), 2) * real_value) {
        if(val > 0) {
          // add it to our set of tagged eigenvalues
          m_tagged_indices.insert(m_tagged_indices.end(), i);
        } else {
          // remove it from the set if it exists
          m_tagged_indices.erase(i);
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_left(const int &val) {
    double real_value = m_shift.real();
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_eigenvalues_alpha.size(); ++i) {
      // if the eigenvalue is to the left of the shift position
      if((m_eigenvalues_alpha[ i ] * std::conj(m_eigenvalues_beta[ i ])).real() <
          std::pow(std::abs(m_eigenvalues_beta[ i ]), 2) * real_value) {
        if(val > 0) {
          // add it to our set of tagged eigenvalues
          m_tagged_indices.insert(m_tagged_indices.end(), i);
        } else {
          // remove it from the set if it exists
          m_tagged_indices.erase(i);
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_upper(const int &val) {
    double imag_value = m_shift.imag();
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_eigenvalues_alpha.size(); ++i) {
      // if the eigenvalue is in the half plane then include it
      if((m_eigenvalues_alpha[ i ] * std::conj(m_eigenvalues_beta[ i ])).imag() >
          std::pow(std::abs(m_eigenvalues_beta[ i ]), 2) * imag_value) {
        if(val > 0) {
          // add it to our set of tagged eigenvalues
          m_tagged_indices.insert(m_tagged_indices.end(), i);
        } else {
          // remove it from the set if it exists
          m_tagged_indices.erase(i);
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_lower(const int &val) {
    double imag_value = m_shift.imag();
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_eigenvalues_alpha.size(); ++i) {
      // if the eigenvalue is in the half plane then include it
      if((m_eigenvalues_alpha[ i ] * std::conj(m_eigenvalues_beta[ i ])).imag() <
          std::pow(std::abs(m_eigenvalues_beta[ i ]), 2) * imag_value) {
        if(val > 0) {
          // add it to our set of tagged eigenvalues
          m_tagged_indices.insert(m_tagged_indices.end(), i);
        } else {
          // remove it from the set if it exists
          m_tagged_indices.erase(i);
        }
      }
    }
  }


  template class DenseLinearEigenSystem<D_complex>
  ;
  template class DenseLinearEigenSystem<double>
  ;

} // end namespace
