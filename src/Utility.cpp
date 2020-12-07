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


  //     template<>
  // double OneD_Node_Mesh<double,double>::maxAbsLocation(unsigned var) {
  //   unsigned N = m_X.size();
  //   return dep_max_abs_location_range(var,m_X[0],m_X[N-1]);
  // }
  
  // template<>
  // double OneD_Node_Mesh<double,double>::dep_max_abs_location_range(unsigned var, double left, double right) {
  //   double max(0.0);
  //   std::size_t maxIndex(0);
  //   // step through the nodes
  //   for(std::size_t node = 0; node < m_X.size(); ++node) {
  //     //std::cout << "m_X[node]=" << m_X[node] << " left=" << left << " right=" << right << "\n";
  //     if ( (m_X[node] >= left) && (m_X[node] <=right) ) {
  //       if(std::abs(m_vars[ node * m_nv + var ]) > max) {
  //         maxIndex = node;
  //         max = std::abs( m_vars[ maxIndex*m_nv + var ]);
  //       }
  //     }
  //   }
  //   if ( ( maxIndex == 0 ) || ( maxIndex == m_X.size()-1 ) ) {
  //     std::cout << "[WARNING] MaxAbsLocationRange: maximumum absolute nodal value is first/last node. \n";
  //     return m_X[ maxIndex ];
  //   }
  //   double f1,f2,f3;
  //   double x1,x2,x3;
  //   f1 = std::abs(m_vars[ (maxIndex-1) * m_nv + var ]);
  //   f2 = std::abs(m_vars[ maxIndex * m_nv + var ]);
  //   f3 = std::abs(m_vars[ (maxIndex+1) * m_nv + var ]);
  //   x1 = m_X[maxIndex-1];
  //   x2 = m_X[maxIndex];
  //   x3 = m_X[maxIndex+1];
  //   return ( f1*(x2+x3)/((x1-x2)*(x1-x3)) + f2*(x1+x3)/((x2-x1)*(x2-x3)) + f3*(x1+x2)/((x3-x1)*(x3-x2)) )
  //     / ( 2.*f1/((x1-x2)*(x1-x3)) + 2.*f2/((x2-x1)*(x2-x3)) + 2.*f3/((x3-x1)*(x3-x2)) );
  // }

    
    double max_abs_location(OneD_Node_Mesh<double> &mesh, unsigned var) {
      double max(0.0);
      std::size_t maxIndex(0);
      // step through the nodes
      for ( std::size_t node = 0; node < mesh.get_nnodes(); ++node ) {
        //std::cout << "m_X[node]=" << m_X[node] << " left=" << left << " right=" << right << "\n";
        if ( std::abs(mesh(node,var)) > max ) {
          maxIndex = node;
          max = std::abs( mesh(node,var) );
        }
      }
      if ( ( maxIndex == 0 ) || ( maxIndex == mesh.get_nnodes()-1 ) ) {
        std::cout << "[WARNING] MaxAbsLocationRange: maximumum absolute nodal value is first/last node. \n";
        return mesh.coord( maxIndex );
      }
      double f1,f2,f3;
      double x1,x2,x3;
      f1 = std::abs(mesh(maxIndex-1,var));
      f2 = std::abs(mesh(maxIndex,var));
      f3 = std::abs(mesh(maxIndex+1,var));
      x1 = mesh.coord(maxIndex-1);
      x2 = mesh.coord(maxIndex);
      x3 = mesh.coord(maxIndex+1);
      return ( f1*(x2+x3)/((x1-x2)*(x1-x3)) + f2*(x1+x3)/((x2-x1)*(x2-x3)) + f3*(x1+x2)/((x3-x1)*(x3-x2)) )
        / ( 2.*f1/((x1-x2)*(x1-x3)) + 2.*f2/((x2-x1)*(x2-x3)) + 2.*f3/((x3-x1)*(x3-x2)) );      
    }


    double max_abs_location_range(OneD_Node_Mesh<double> &mesh, unsigned var, double left, double right ) {
      double max(0.0);
      std::size_t maxIndex(0);
      // step through the nodes
      for ( std::size_t node = 0; node < mesh.get_nnodes(); ++node ) {
        //std::cout << "m_X[node]=" << m_X[node] << " left=" << left << " right=" << right << "\n";
        if ( ( mesh.coord(node) >= left ) && ( mesh.coord(node) <= right ) ) {
          if ( std::abs(mesh(node,var)) > max ) {
            maxIndex = node;
            max = std::abs( mesh(node,var) );
          }
        }
      }
      if ( ( maxIndex == 0 ) || ( maxIndex == mesh.get_nnodes()-1 ) ) {
        std::cout << "[WARNING] MaxAbsLocationRange: maximumum absolute nodal value is first/last node. \n";
        return mesh.coord( maxIndex );
      }
      double f1,f2,f3;
      double x1,x2,x3;
      f1 = std::abs(mesh(maxIndex-1,var));
      f2 = std::abs(mesh(maxIndex,var));
      f3 = std::abs(mesh(maxIndex+1,var));
      x1 = mesh.coord(maxIndex-1);
      x2 = mesh.coord(maxIndex);
      x3 = mesh.coord(maxIndex+1);
      return ( f1*(x2+x3)/((x1-x2)*(x1-x3)) + f2*(x1+x3)/((x2-x1)*(x2-x3)) + f3*(x1+x2)/((x3-x1)*(x3-x2)) )
        / ( 2.*f1/((x1-x2)*(x1-x3)) + 2.*f2/((x2-x1)*(x2-x3)) + 2.*f3/((x3-x1)*(x3-x2)) );      
    }

    

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
