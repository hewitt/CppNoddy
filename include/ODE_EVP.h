/// \file ODE_EVP.h
/// A specification of a class for an \f$ n^{th} \f$-order ODE LINEAR EVP defined by
/// \f[ \lambda M_1({\underline f}(x),x ) \cdot {\underline f } + M_0({\underline f}(x),x ) \cdot {\underline f}^\prime (x) = {\underline R}( {\underline f}(x), x )\,, \f]
/// subject to \f$ n \f$ zero Dirichlet conditions defined at \f$ x = x_{left} \f$ or
/// \f$ x_{right} \f$ for some components of \f$ {\underline f}(x) \f$. The routine
/// constructs a banded 2nd order finite-difference (banded) representation of the EVP
/// in the form \f[ A_{M\times M} {\underline x} = \lambda B_{M\times M} {\underline x} \f]
/// which can then be solved via the LinearEigenSystem class where
/// \f$ M = N \times O \f$ for an order O equation and N (uniformly spaced) mesh points. The default
/// configuration uses the DenseLinearEigenSystem class for solution via LAPACK
/// although ARPACK can be used if you really know the resulting system has the
/// correct mass matrix properties (typically it wont!).

#ifndef ODE_EVP_H
#define ODE_EVP_H

#include <memory>
#include <algorithm>

#include <DenseVector.h>
#include <DenseMatrix.h>
#include <Equation_2matrix.h>
#include <Exceptions.h>
#include <Uncopyable.h>
#include <OneD_Node_Mesh.h>
#include <LinearEigenSystem_base.h>
#include <Residual.h>

namespace CppNoddy
{

  /// A templated object for real/complex vector system
  /// of first-order ordinary differential equations.

  template <typename _Type>
  class ODE_EVP : private Uncopyable
  {

  public:

    /// The class is defined by a vector function for the system.
    /// \param equation_ptr A pointer to an equation with 2 associated matrices; matrix1 will define the eigenvalue problem.
    /// \param nodes A vector of nodal points.
    /// \param ptr_to_left_residual A pointer to a residual object that defines the LHS boundary conditions.
    /// \param ptr_to_right_residual A pointer to a residual object that defines the RHS boundary conditions.
    ODE_EVP( Equation_2matrix<_Type > *equation_ptr,
             const DenseVector<double> &nodes,
             Residual<_Type>* ptr_to_left_residual,
             Residual<_Type>* ptr_to_right_residual );

    /// Destructor
    ~ODE_EVP();

    /// Formulate and solve the global eigenvalue problem
    /// for a linear system.
    void eigensolve();

    /// Allow access to the underlying dense linear eigensystem
    /// through a pointer to the private member data.
    LinearEigenSystem_base *p_eigensystem();

    void add_tagged_to_mesh()
    {
      // clear the existing data if any
      MESHES.clear();
      // order of the equation
      unsigned order = p_EQUATION -> get_order();
      // get the eigenvalues
      DenseVector<D_complex> vals( p_SYSTEM -> get_tagged_eigenvalues() );
      // get the eigenvectors
      DenseMatrix<D_complex> vecs( p_SYSTEM -> get_tagged_eigenvectors() );
      // loop through the eigenvectors
      for ( unsigned ivec = 0; ivec < vals.size(); ++ivec )
      {
        // make a mesh with the right node distribution
        // we'll increase the order by 1 to allow the eigenvalue to be
        // stored at each nodal point -- this is wasteful but very useful
        // for feeding this into a BVP for local refinement.
        OneD_Node_Mesh<D_complex> eigfn( NODES, order + 1 );
        // loop through all nodes
        for ( unsigned node = 0; node < NODES.size(); ++node )
        {
          // complex vector of the dof at this node ( + 1 for the eigenvalue)
          DenseVector<D_complex> vars_at_node( order + 1, 0.0 );
          // get the dof from the eigenvector
          for ( unsigned var = 0; var < order; ++var )
          {
            vars_at_node[ var ] = vecs[ ivec ][ node * order + var ];
          }
          // the last variable at each node is the corresponding eigenvalue
          vars_at_node[ order ] = vals[ ivec ];
          // set the first 'order' dof to be those from the eigenvector
          eigfn.set_nodes_vars( node, vars_at_node );
        }
        //// store the eigenvalue in the mesh at each node ... wasteful, but useful
        //// a complex vector filled with the same value N times
        //DenseVector<D_complex> c( NODES.size(), vals[ ivec ] );
        //// add it to the mesh -- for use in nonlinear local refinement via ODE_BVP
        //eigfn.push_var( c );
        // add the eigenfunction to the vector of meshes
        MESHES.push_back( eigfn );
      }
    }

    OneD_Node_Mesh<D_complex> get_mesh( const unsigned& i ) const
    {
#ifdef PARANOID
      if ( i > MESHES.size() )
      {
        std::string problem;
        problem = "You have tried to extract an eigenfunction from the ODE_EVP class\n";
        problem += "whose index is outside the range of stored meshes.\n";
        throw ExceptionRange( problem, MESHES.size(), i );
      }
#endif
      return MESHES[ i ];
    }

  private:

    /// A method that constructs the banded matrix problem
    void assemble_dense_problem();

    /// The function associated with this instance.
    Equation_2matrix<_Type > *p_EQUATION;
    /// Pointer to the residual defining the LHS BC
    Residual<_Type > *p_LEFT_RESIDUAL;
    /// Pointer to the residual defining the RHS BC
    Residual<_Type > *p_RIGHT_RESIDUAL;
    /// The dense linear eigensystem
    LinearEigenSystem_base *p_SYSTEM;
    /// Matrices for the LAPACK QZ routine -- must be DENSE
    DenseMatrix<_Type>* p_A_DENSE;
    DenseMatrix<_Type>* p_B_DENSE;
    /// the nodal distribution
    DenseVector<double> NODES;
    /// A vector of uniform meshes that store the eigenfunctions
    /// and eigenvalue (as the last dof at each nodal point)
    /// for use in later local refinement via the ODE_BVP class
    std::vector< OneD_Node_Mesh<D_complex> > MESHES;
    /// has the eigensystem been constructed/solved?
    bool CONSTRUCTED;

  }
  ; // end class


} // end namespace

#endif

