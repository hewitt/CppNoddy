/// \file LinearEigenSystem_base.cpp
/// Implementation for the LinearEigenSystem_base class.
/// The specific eigesolvers inherit from here.

#include <vector>
#include <set>

#include <LinearEigenSystem_base.h>
#include <Exceptions.h>
#include <Types.h>

namespace CppNoddy {

  LinearEigenSystem_base::LinearEigenSystem_base() :
    m_shift(D_complex(0., 0.)),
    m_calc_eigenvectors(true) {
  }


  LinearEigenSystem_base::~LinearEigenSystem_base()
  {}


  void LinearEigenSystem_base::set_shift(const D_complex& z) {
    m_shift = z;
  }


  std::complex<double> LinearEigenSystem_base::get_shift() const {
    return m_shift;
  }


  void LinearEigenSystem_base::eigensolve() {
    std::string problem;
    problem = "The LinearEigenSystem_base::eigensolve method has been called\n";
    problem += "but the method has not been implemented ... this should be \n";
    problem += "implemented in the sub-class.";
    throw ExceptionExternal(problem);
  }


  void LinearEigenSystem_base::set_calc_eigenvectors(bool flag) {
    m_calc_eigenvectors = flag;
  }


  DenseVector<D_complex> LinearEigenSystem_base::get_tagged_eigenvalues() const {
    if(m_tagged_indices.size() == 0) {
      std::string problem;
      problem = "In LinearEigenSystem_base.get_tagged_eigenvalues() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime(problem);
    }
    // storage for the eigenvalues
    DenseVector<D_complex> evals;
    // loop through the tagged set
    for(iter p = m_tagged_indices.begin(); p != m_tagged_indices.end(); ++p) {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      std::cout << " number " << j << " is a tagged ev.\n";
      // work out the complex eigenvalue associated with this index
      // and add it to the vector
      evals.push_back(m_all_eigenvalues[ j ]);
    }
    // return the complex vector of eigenvalues
    return evals;
  }


  DenseMatrix<D_complex> LinearEigenSystem_base::get_tagged_eigenvectors() const {
    if(m_tagged_indices.size() == 0) {
      std::string problem;
      problem = "In LinearEigenSystem_base.get_tagged_eigenvectors() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime(problem);
    }
    // number of degrees of freedom in each eigenvector
    std::size_t n = m_all_eigenvectors[0].size();
    // eigenvector storage : size() eigenvectors each of length n
    DenseMatrix<D_complex> evecs(m_tagged_indices.size(), n, 0.0);
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


  void LinearEigenSystem_base::tag_eigenvalues_all() {
    for(unsigned i = 0; i < m_all_eigenvalues.size(); ++i) {
      m_tagged_indices.insert(m_tagged_indices.end(), i);
    }
    std::cout << "** I have tagged " << m_all_eigenvalues.size() << " eigenvalues.\n";
  }


  void LinearEigenSystem_base::tag_eigenvalues_disc(const int &val, const double& radius) {
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_all_eigenvalues.size(); ++i) {
      // if the eigenvalue is in the disc centred at shift then include it
      if(std::abs(m_all_eigenvalues[ i ] - m_shift) < radius) {
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


  void LinearEigenSystem_base::tag_eigenvalues_right(const int &val) {
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_all_eigenvalues.size(); ++i) {
      // if the eigenvalue is in the disc centred at m_shift then include it
      if((m_all_eigenvalues[ i ] - m_shift).real() > 0.0) {
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


  void LinearEigenSystem_base::tag_eigenvalues_left(const int &val) {
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_all_eigenvalues.size(); ++i) {
      // if the eigenvalue is in the disc centred at m_shift then include it
      if((m_all_eigenvalues[ i ] - m_shift).real() < 0.0) {
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


  void LinearEigenSystem_base::tag_eigenvalues_upper(const int &val) {
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_all_eigenvalues.size(); ++i) {
      // if the eigenvalue is in the disc centred at m_shift then include it
      if((m_all_eigenvalues[ i ] - m_shift).imag() > 0.0) {
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

  void LinearEigenSystem_base::tag_eigenvalues_lower(const int &val) {
    // loop through all the eigenvalues
    for(std::size_t i = 0; i < m_all_eigenvalues.size(); ++i) {
      // if the eigenvalue is in the disc centred at m_shift then include it
      if((m_all_eigenvalues[ i ] - m_shift).imag() < 0.0) {
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


  class LinearEigenSystem_base;

} // end namespace
