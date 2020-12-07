/// \file Utility.h
/// A spec for a collection of utility functions.

#ifndef UTILITY_H
#define UTILITY_H

#include <OneD_Node_Mesh.h>
#include <string>
#include <ctime>
#include <numeric>
#include <algorithm>
#include <numeric>
#include <sys/types.h>
#include <sys/stat.h>

#include <Types.h>
#include <Exceptions.h>
#include <FortranBLAS.h>
#include <TwoD_Node_Mesh.h>
#include <TwoD_Mapped_Node_Mesh.h>
#include <Sequential_Matrix_base.h>

namespace CppNoddy {
  /// Some utility methods associated with CppNoddy containers.
  namespace Utility {

    /// Return a DENSE vector with the nodal points of a uniform
    /// mesh distributed between the upper/lower bounds as specified
    /// \param lower The lower bound of the uniform nodal distribution
    /// \param upper The upper bound of the uniform nodal distribution
    /// \param N The number of nodal points
    DenseVector<double> uniform_node_vector(const double& lower, const double& upper, const std::size_t& N);

    /// Return a DENSE vector with the nodal points of a non-uniform
    /// mesh distributed between the upper/lower bounds as specified
    /// with more nodes clustered near lower or upper depending upon
    /// the differencee of the power from unity. When power=1 this should
    /// provide a uniform mesh.
    /// \param lower The lower bound of the uniform nodal distribution
    /// \param upper The upper bound of the uniform nodal distribution
    /// \param N The number of nodal points
    /// \param power A measure of the non-uniformity
    /// \return A vector of nodal positions with a power law distribution
    DenseVector<double> power_node_vector(const double& lower, const double& upper, const std::size_t& N, const double& power);

    /// Return a dense vector with two uniform distributions in two separate
    /// regions.
    /// \param lower The first node
    /// \param mid The node that defines the boundary between the uniform meshes
    /// \param upper The final node
    /// \param N1 The number of nodes in the first region
    /// \param N2 The number of nodes in the second region
    /// \return A combined vector of nodes of length N1+N2
    DenseVector<double> two_uniform_node_vector(const double& lower, const double& mid, const double& upper, const std::size_t& N1, const std::size_t& N2);

    /// Return a dense vector with two uniform distributions in two separate
    /// regions.
    /// \param lower The first node
    /// \param mid1 The node that defines the first interior boundary
    /// \param mid2 The node that defines the second interior boundary
    /// \param upper The final node
    /// \param N1 The number of nodes in the first region
    /// \param N2 The number of nodes in the second region
    /// \param N3 The number of nodes in the third region
    /// \return A combined vector of nodes of length N1+N2+N3
    DenseVector<double> three_uniform_node_vector(const double& lower, const double& mid1, const double& mid2, const double& upper, const std::size_t& N1, const std::size_t& N2, const std::size_t& N3);

    /// Return a dense vector of nodal positions with more nodes concentrated
    /// at the mid point of the range.
    /// \param lower The first nodal position.
    /// \param upper The final nodal position.
    /// \param N The number of nodes required.
    /// \param power A measure of the non-uniformity, power = 1 => uniform distribution
    DenseVector<double> mid_weighted_node_vector(const double& lower, const double& upper, const std::size_t& N, const double& power);

    //
    //
    // SOME typical ops
    //
    //
    
    template <typename _Type>
    void vels_from_streamfn_Cartesian(const TwoD_Node_Mesh<_Type>& source, TwoD_Node_Mesh<_Type>& uv) {
      std::cout << "PLAIN MESH version\n";
      std::size_t Nx(source.get_nnodes().first);
      std::size_t Ny(source.get_nnodes().second);
      double dx(source.coord(1,1).first - source.coord(0,0).first);
      double dy(source.coord(1,1).second - source.coord(0,0).second);
      // differentiate the streamfunction to get the velocity field
      {
        // west internal nodes
        std::size_t i(0);
        for(std::size_t j = 1; j < Ny - 1; ++j) {
          uv(i, j, 0) = (source(i, j + 1, 0) - source(i, j - 1, 0)) / (2 * dy);
          uv(i, j, 1) = -(-source(i + 2, j, 0) + 4. * source(i + 1, j, 0) - 3. * source(i, j, 0)) / (2 * dx);
        }
      }
      {
        // east internal nodes
        std::size_t i(Nx - 1);
        for(std::size_t j = 1; j < Ny - 1; ++j) {
          uv(i, j, 0) = (source(i, j + 1, 0) - source(i, j - 1, 0)) / (2 * dy);
          uv(i, j, 1) = -(source(i - 2, j, 0) - 4. * source(i - 1, j, 0) + 3. * source(i, j, 0)) / (2 * dx);
        }
      }
      {
        // south internal nodes
        std::size_t j(0);
        for(std::size_t i = 1; i < Nx - 1; ++i) {
          uv(i, j, 0) = (-source(i, j + 2, 0) + 4. *  source(i, j + 1, 0) - 3. * source(i, j, 0)) / (2 * dy);
          uv(i, j, 1) = -(source(i + 1, j, 0) - source(i - 1, j, 0)) / (2 * dx);
        }
      }
      {
        // north internal nodes
        std::size_t j(Ny - 1);
        for(std::size_t i = 1; i < Nx - 1; ++i) {
          uv(i, j, 0) = (source(i, j - 2, 0) - 4. *  source(i, j - 1, 0) + 3. * source(i, j, 0)) / (2 * dy);
          uv(i, j, 1) = -(source(i + 1, j, 0) - source(i - 1, j, 0)) / (2 * dx);
        }
      }
      {
        // corner nodes
        {
          // sw
          std::size_t i(0);
          std::size_t j(0);
          uv(i, j, 0) = (-source(i, j + 2, 0) + 4. *  source(i, j + 1, 0) - 3. * source(i, j, 0)) / (2 * dy);
          uv(i, j, 1) = -(-source(i + 2, j, 0) + 4. * source(i + 1, j, 0) - 3. * source(i, j, 0)) / (2 * dx);
        }
        {
          // nw
          std::size_t i(0);
          std::size_t j(Ny - 1);
          uv(i, j, 0) = (source(i, j - 2, 0) - 4. *  source(i, j - 1, 0) + 3. * source(i, j, 0)) / (2 * dy);
          uv(i, j, 1) = -(-source(i + 2, j, 0) + 4. * source(i + 1, j, 0) - 3. * source(i, j, 0)) / (2 * dx);
        }
        {
          // ne
          std::size_t i(Nx - 1);
          std::size_t j(Ny - 1);
          uv(i, j, 0) = (source(i, j - 2, 0) - 4. *  source(i, j - 1, 0) + 3. * source(i, j, 0)) / (2 * dy);
          uv(i, j, 1) = -(source(i - 2, j, 0) - 4. * source(i - 1, j, 0) + 3. * source(i, j, 0)) / (2 * dx);
        }
        {
          // se
          std::size_t i(Nx - 1);
          std::size_t j(0);
          uv(i, j, 0) = (-source(i, j + 2, 0) + 4. *  source(i, j + 1, 0) - 3. * source(i, j, 0)) / (2 * dy);
          uv(i, j, 1) = -(source(i - 2, j, 0) - 4. * source(i - 1, j, 0) + 3. * source(i, j, 0)) / (2 * dx);
        }
      }
      {
        // interior nodes
        for(std::size_t i = 1; i < Nx - 1; ++i) {
          for(std::size_t j = 1; j < Ny - 1; ++j) {
            uv(i, j, 0) = (source(i, j + 1, 0) - source(i, j - 1, 0)) / (2 * dy);
            uv(i, j, 1) = -(source(i + 1, j, 0) - source(i - 1, j, 0)) / (2 * dx);
          }
        }
      }
    }


    template <typename _Type>
    void vels_from_streamfn_Cartesian(const TwoD_Mapped_Node_Mesh<_Type>& source, TwoD_Mapped_Node_Mesh<_Type>& uv) {
      std::cout << "MAPPED MESH version\n";
      std::size_t Nx(source.get_nnodes().first);
      std::size_t Ny(source.get_nnodes().second);
      double dX( source.get_comp_step_sizes().first );
      double dY( source.get_comp_step_sizes().second );
      // differentiate the streamfunction to get the velocity field
      {
        // west internal nodes
        std::size_t i(0);
        for(std::size_t j = 1; j < Ny - 1; ++j) {
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);
          uv(i, j, 0) = (source(i, j + 1, 0) - source(i, j - 1, 0)) * Yd / (2 * dY);
          uv(i, j, 1) = -(-source(i + 2, j, 0) + 4. * source(i + 1, j, 0) - 3. * source(i, j, 0)) * Xd / (2 * dX);
        }
      }
      {
        // east internal nodes
        std::size_t i(Nx - 1);
        for(std::size_t j = 1; j < Ny - 1; ++j) {
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);
          uv(i, j, 0) = (source(i, j + 1, 0) - source(i, j - 1, 0)) * Yd / (2 * dY);
          uv(i, j, 1) = -(source(i - 2, j, 0) - 4. * source(i - 1, j, 0) + 3. * source(i, j, 0)) * Xd / (2 * dX);
        }
      }
      {
        // south internal nodes
        std::size_t j(0);
        for(std::size_t i = 1; i < Nx - 1; ++i) {
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);
          uv(i, j, 0) = (-source(i, j + 2, 0) + 4. *  source(i, j + 1, 0) - 3. * source(i, j, 0))* Yd / (2 * dY);
          uv(i, j, 1) = -(source(i + 1, j, 0) - source(i - 1, j, 0)) * Xd / (2 * dX);
        }
      }
      {
        // north internal nodes
        std::size_t j(Ny - 1);
        for(std::size_t i = 1; i < Nx - 1; ++i) {
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);
          uv(i, j, 0) = (source(i, j - 2, 0) - 4. *  source(i, j - 1, 0) + 3. * source(i, j, 0)) * Yd / (2 * dY);
          uv(i, j, 1) = -(source(i + 1, j, 0) - source(i - 1, j, 0)) * Xd / (2 * dX);
        }
      }
      {
        // corner nodes
        {
          // sw
          std::size_t i(0);
          std::size_t j(0);
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);	  
          uv(i, j, 0) = (-source(i, j + 2, 0) + 4. *  source(i, j + 1, 0) - 3. * source(i, j, 0))* Xd / (2 * dY);
          uv(i, j, 1) = -(-source(i + 2, j, 0) + 4. * source(i + 1, j, 0) - 3. * source(i, j, 0))* Yd / (2 * dY);
        }
        {
          // nw
          std::size_t i(0);
          std::size_t j(Ny - 1);
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);	  
          uv(i, j, 0) = (source(i, j - 2, 0) - 4. *  source(i, j - 1, 0) + 3. * source(i, j, 0)) * Yd / (2 * dY);
          uv(i, j, 1) = -(-source(i + 2, j, 0) + 4. * source(i + 1, j, 0) - 3. * source(i, j, 0))* Xd / (2 * dX);
        }
        {
          // ne
          std::size_t i(Nx - 1);
          std::size_t j(Ny - 1);
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);	  
          uv(i, j, 0) = (source(i, j - 2, 0) - 4. *  source(i, j - 1, 0) + 3. * source(i, j, 0)) * Yd / (2 * dY);
          uv(i, j, 1) = -(source(i - 2, j, 0) - 4. * source(i - 1, j, 0) + 3. * source(i, j, 0)) * Xd/ (2 * dX);
        }
        {
          // se
          std::size_t i(Nx - 1);
          std::size_t j(0);
	  double x = source.coord(i,j).first;
	  double Xd = source.FnComp_Xd(x);
	  double y = source.coord(i,j).second;
	  double Yd = source.FnComp_Yd(y);	  
          uv(i, j, 0) = (-source(i, j + 2, 0) + 4. *  source(i, j + 1, 0) - 3. * source(i, j, 0))* Yd / (2 * dY);
          uv(i, j, 1) = -(source(i - 2, j, 0) - 4. * source(i - 1, j, 0) + 3. * source(i, j, 0)) * Xd / (2 * dX);
        }
      }
      {
        // interior nodes
        for(std::size_t i = 1; i < Nx - 1; ++i) {
          for(std::size_t j = 1; j < Ny - 1; ++j) {
	    double x = source.coord(i,j).first;
	    double Xd = source.FnComp_Xd(x);
	    double y = source.coord(i,j).second;
	    double Yd = source.FnComp_Yd(y);	  
            uv(i, j, 0) = (source(i, j + 1, 0) - source(i, j - 1, 0)) * Yd / (2 * dY);
            uv(i, j, 1) = -(source(i + 1, j, 0) - source(i - 1, j, 0)) * Xd / (2 * dX);
          }
        }
      }
    }



    
    //
    //
    // MATRIX OPERATIONS
    //
    //

    /// BLAS wrapper to do DOUBLE DENSE A_{MxK} * B_{KxN} = C_{MxN}
    /// Since this is a Fortran library, it assumes a column_major
    /// format, but CppNoddy uses row_major. To be consistent we'll
    /// simply do (B^T)_{NxK} * (A^T)_{KxM} = (C^T)_{NxM} instead.
    /// Note that inversion of the transpose of the result C^T
    /// is handled implicitly via the construction of C.
    /// \param A First dense double matrix to be multiplied
    /// \param B Second dense double matrix to be multiplied
    /// \return The result of the multiplication C=A*B
    DenseMatrix<double> multiply(DenseMatrix<double>& A, DenseMatrix<double>& B);

    //
    //
    // VECTOR OPERATIONS
    //
    //

    /// Templated dot product.
    /// \param X First dense vector
    /// \param Y Second dense vector
    /// \return The dot product
    template <typename _Type>
    _Type dot(const DenseVector<_Type>& X, const DenseVector<_Type>& Y) {
      if(X.size() != Y.size()) {
        std::string problem;
        problem = "The Utilities::dot method has been called \n";
        problem += "with two unequal length vectors.";
        throw ExceptionGeom(problem, X.size(), Y.size());
      }
      return inner_product(X.begin(), X.end(), Y.begin(), _Type(0.0));
    }

    template <typename _Type>
    int sgn(const _Type& a) {
      if(a > (_Type)0) {
        return 1;
      } else if(a < (_Type)0) {
        return -1;
      } else {
        return 0;
      }
    }

    //
    //
    // COMPLEX UTILS
    //
    //

    /// Return a double DENSE vector containing the real part
    /// of a complex DENSE vector
    /// \param X The complex vector to take the real part of
    DenseVector<double> real(const DenseVector<D_complex>& X);

    /// Return a double DENSE vector containing the imaginary part
    /// of a complex DENSE vector
    /// \param X The complex vector to take the imaginary part of
    DenseVector<double> imag(const DenseVector<D_complex>& X);

    //
    //
    // STRING TWEAKERY
    //
    //

    /// Return an integer value as a string - useful for file naming
    /// \param val The integer value to be stringified.
    std::string stringify(const int &val);

    /// Return a double value as a string - useful for file naming.
    /// \param val The double value to be stringified
    /// \param p Precision to be used in the output
    std::string stringify(const double &val, int p);


    
    double max_abs_location( OneD_Node_Mesh<double> &mesh , unsigned var);
    
    double max_abs_location_range( OneD_Node_Mesh<double> &mesh , unsigned var, double left, double right);
    
  }

} // end namespace

#endif // UTILITY_H
