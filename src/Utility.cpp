/// \file Utility.cpp
/// An implementation for a collection of utility functions.

#include <string>
#include <ctime>

#include <Types.h>
#include <Exceptions.h>
#include <Utility.h>
#include <FortranData.h>
#include <FortranBLAS.h>

namespace CppNoddy {
  namespace Utility {
    
    DenseVector<double> uniform_node_vector(const double& lower, const double& upper, const std::size_t& N) {
      DenseVector<double> V;
      V.reserve(N);
      const double delta = (upper - lower) / (N - 1);
      for(std::size_t i = 0; i < N; ++i) {
        V.push_back(lower + delta * i);
      }
      return V;
    }

    DenseVector<double> power_node_vector(const double& lower, const double& upper, const std::size_t& N, const double& power) {
      DenseVector<double> V(N, 0.0);
      for(std::size_t i = 0; i < N; ++i) {
        V[ i ] = lower + (upper - lower) * std::pow((double)i / (N - 1), power);
      }
      return V;
    }

    DenseVector<double> two_uniform_node_vector(const double& lower, const double& mid, const double& upper, const std::size_t& N1, const std::size_t& N2) {
      DenseVector<double> first = uniform_node_vector(lower, mid, N1);
      DenseVector<double> second = uniform_node_vector(mid, upper, N2 + 1);
      // skip the common elt by starting at 1 not 0
      for(std::size_t i = 1; i < N2+1; ++i) {
        first.push_back(second[ i ]);
      }
      return first;
    }

    DenseVector<double> three_uniform_node_vector(const double& lower, const double& mid1, const double& mid2, const double& upper, const std::size_t& N1, const std::size_t& N2, const std::size_t& N3) {
      DenseVector<double> first = uniform_node_vector(lower, mid1, N1);
      DenseVector<double> second = uniform_node_vector(mid1, mid2, N2+1);
      DenseVector<double> third = uniform_node_vector(mid2, upper, N3+1);
      // skip the common elt by starting at 1 not 0
      for(std::size_t i = 1; i < N2+1; ++i) {
        first.push_back(second[ i ]);
      }
      for(std::size_t i = 1; i < N3+1; ++i) {
        first.push_back(third[ i ]);
      }
      return first;
    }

    DenseVector<double> mid_weighted_node_vector(const double& lower, const double& upper, const std::size_t& N, const double& weight) {
      DenseVector<double> node_vector(N, 0.0);
      // make a center weighted distribution over -1 to 1
      for(std::size_t i = 0; i < N; ++i) {
        double s(-1.0 + 2.0 * i / (N - 1));
        node_vector[ i ] = (weight * std::pow(s, 3) + s) / (weight + 1);
      }
      // map the -1 to 1 range to lower to upper
      for(std::size_t i = 0; i < N; ++i) {
        // move the range to 0->2
        node_vector[ i ] += 1.0;
        // move the range tp 0->1
        node_vector[ i ] /= 2.0;
        // move the range to lower -> upper
        node_vector[ i ] *= (upper - lower);
        node_vector[ i ] += lower;
      }
      return node_vector;
    }



    DenseVector<double> real(const DenseVector<D_complex>& X) {
      DenseVector<double> temp(X.size(), 0.0);
      for(std::size_t i = 0; i < X.size(); ++i) {
        temp[ i ] = X[ i ].real();
      }
      return temp;
    }

    DenseVector<double> imag(const DenseVector<D_complex>& X) {
      DenseVector<double> temp(X.size(), 0.0);
      for(std::size_t i = 0; i < X.size(); ++i) {
        temp[ i ] = X[ i ].imag();
      }
      return temp;
    }

    std::string stringify(const int& val) {
      std::stringstream temp;
      temp << val;
      return temp.str();
    }

    std::string stringify(const double& val, int p) {
      std::stringstream temp;
      temp.precision(p);
      temp << val;
      return temp.str();
    }

    //        // dot product specialised for doubles to use
    //        // the BLAS library if LAPACK specified
    //        double dot<double>( const DenseVector<double>& X, const DenseVector<double>& Y )
    //        {
    //            // check lengths
    //            if ( X.size() != Y.size() )
    //            {
    //                std::string problem( " The LAPACK::dot method has detected a geometry error. \n" );
    //                throw ExceptionGeom( problem, X.size(), Y.size() );
    //            }
    //#ifndef LAPACK
    //            double sum( 0.0 );
    //            inner_product( X.begin(), X.end(), Y.begin(), sum );
    //            return sum;
    //#else
    //            return BLAS_DDOT( X.size(), &X[0], 1, &Y[0], 1 );
    //#endif
    //        }

    DenseMatrix<double> multiply(DenseMatrix<double>& A, DenseMatrix<double>& B) {
#ifndef LAPACK
      std::string problem;
      problem = "The Utilities::multiply method has been called\n";
      problem += "but the compiler option -DLAPACK was not provided when\n";
      problem += "the library was built. This non-member function requires BLAS.";
      throw ExceptionExternal(problem);
#else
      // set the matrix geometries
      std::size_t M = A.nrows();
      std::size_t N = B.ncols();
      std::size_t K = A.ncols();
      // No need to transpose first, because we will in fact do B^T * A^T = C^T
      // New the memory for the result
      FortranData Af(A, false);
      FortranData Bf(B, false);
      FortranData Cf(M * N);
#ifdef PARANOID

      if(K != B.nrows()) {
        std::string problem(" The LAPACK::multiply method has detected a failure \n");
        throw ExceptionGeom(problem, A.nrows(), A.ncols(), B.nrows(), B.ncols());
      }
#endif
      // call Fortran BLAS
      // call Fortran BLAS
      BLAS_DGEMM((char*) "N", (char*) "N", N, M, K, 1.0, Bf.base(), N, Af.base(), K, 0.0, Cf.base(), N);
      // Return a DenseMatrix<double> from the results -- since Cf is in column_major
      // format, this will actually transpose the Cf data to provide C as required
      return (Cf.to_dense_matrix(M, N));
#endif

    }

  }

} // end namespace
