/// \file ODE_EVP.cpp
/// A specification of a class for an \f$ n^{th} \f$-order ODE LINEAR EVP defined by
/// \f[ {\underline f}^\prime (x) = {\underline R}( {\underline f}(x), x )\,, \f]
/// subject to \f$ n \f$ zero Dirichlet conditions defined at \f$ x = x_{left} \f$ or
/// \f$ x_{right} \f$ for some components of \f$ {\underline f}(x) \f$. The routine
/// constructs a banded 2nd order finite-difference (banded) representation of the EVP
/// in the form \f[ A_{MxM} {\underline x} = \lambda B_{MxM} {\underline x} \f]
/// which can then be solved via the LinearEigenSystem class where
/// \f$ M = N x O \f$ for an order O equation and N (uniformly spaced) mesh points. The default
/// configuration uses the DenseLinearEigenSystem class for solution via LAPACK
/// although ARPACK can be used if you really know the resulting system has the
/// correct mass matrix properties (typically it wont!).

#include <string>

#include <Types.h>
#include <Equation.h>
#include <Utility.h>
#include <ODE_EVP.h>
#include <Exceptions.h>
#include <DenseLinearEigenSystem.h>
#include <BandedLinearEigenSystem.h>

#include <Equation_with_mass.h>

namespace CppNoddy
{
  template <typename _Type>
  ODE_EVP<_Type>::ODE_EVP( Equation_with_mass<_Type >* ptr_to_equation,
                           const DenseVector<double> &nodes,
                           Residual<_Type>* ptr_to_left_residual,
                           Residual<_Type>* ptr_to_right_residual,
                           std::string which ) :
      p_EQUATION( ptr_to_equation ),
      p_LEFT_RESIDUAL( ptr_to_left_residual ),
      p_RIGHT_RESIDUAL( ptr_to_right_residual ),
      NODES( nodes ),  
      VERSION( which )
  {
    unsigned n( nodes.size() );      
    unsigned order( p_EQUATION -> get_order() );
    // banded storage of the eigenproblem
    p_A = new BandedMatrix<_Type>( n * order, 2 * order - 1, 0.0 );
    p_B = new BandedMatrix<_Type>( n * order, 2 * order - 1, 0.0 );
   // construct the banded matrix eigenvalue problem that corresponds
    // to the equation & boundary conditions specified in the constructor.
    construct_banded_problem();
    // now we have the problem, we just need to pick a method of solution
    if ( "lapack" == VERSION )
    {
      // user wants to solve the dense evp via LAPACK QZ routine
      // first make the dense matrix equivalent of the problem
      // using the Dense( banded ) constructor
      p_A_DENSE = new DenseMatrix<_Type>( *p_A );
      p_B_DENSE = new DenseMatrix<_Type>( *p_B );
      // point the system pointer to the dense problem
      p_SYSTEM = new DenseLinearEigenSystem<_Type>( p_A_DENSE, p_B_DENSE );
    }
    else
    {
      if ( "arpack" == VERSION )
      {
        // ARPACK solve is easy, since the problem is already
        // in the appropriate banded form
        p_SYSTEM = new BandedLinearEigenSystem<_Type>( p_A, p_B );
        std::cout << " [Warning] : using ARPACK as the solver in the ODE_EVP routine\n";
        std::cout << " is probably not a good idea. The formulation of the matrix \n";
        std::cout << " has to be done in such a way that the mass matrix is symmetric \n";
        std::cout << " positive semi-definite. Please use the slower LAPACK QZ algorithm.\n";
      }
      else
      {
        std::string problem;
        problem = "The ODE_EVP::eigensolve method has been called with a request for \n";
        problem += "an eigensolver method that has not been implemented. Currently\n";
        problem += "you must use either 'lapack' or 'arpack'\n";
        throw ExceptionRuntime( problem );
      }
    }     
 }


  template <typename _Type>
  ODE_EVP<_Type>::~ODE_EVP()
  {
    // clean up the matrices & evp system
    delete p_A;
    delete p_B;
    delete p_SYSTEM;
    // if we used lapack we also have some dense matrices to clean up
    if ( "lapack" == VERSION )
    {
      delete p_A_DENSE;
      delete p_B_DENSE;
    }
  }

  template <typename _Type>
  LinearEigenSystem_base* ODE_EVP<_Type>::p_eigensystem()
  {
    return p_SYSTEM;
  }

  template <typename _Type>
  void ODE_EVP<_Type>::eigensolve()
  {
    // solve the system
    p_SYSTEM -> eigensolve();
  }


  template <typename _Type>
  void ODE_EVP<_Type>::construct_banded_problem()
  {
    unsigned order( p_EQUATION -> get_order() );
    // Jacobian matrix for the equation
    DenseMatrix<_Type> jac_midpt( order, order, 0.0 );
    // local state variable and functions
    // the problem has to be linear, so this is a dummy vector.
    DenseVector<_Type> temp_dummy( order, 0.0 );
    // number of nodes in the mesh      
    std::size_t num_of_nodes( NODES.size() );      
    // row counter
    std::size_t row( 0 );

    // update the BC residuals for the current iteration
    p_LEFT_RESIDUAL -> update( temp_dummy );
    // add the (linearised) LHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_LEFT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        ( *p_A ) ( row, var ) = p_LEFT_RESIDUAL -> jacobian()( i, var );
        ( *p_B ) ( row, var ) = 0.0;
      }
      ++row;
    }

    // inner nodes of the mesh, node = 0,1,2,...,num_of_nodes-2
    for ( std::size_t node = 0; node <= num_of_nodes - 2; ++node )
    {
      const std::size_t lnode = node;
      const std::size_t rnode = node + 1;
      // get the current step length
      double h = NODES[ rnode ] - NODES[ lnode ];
      // reciprocal of the spatial step
      const double invh( 1. / h );
      // mid point of the independent variable
      double x_midpt = 0.5 * ( NODES[ lnode ] + NODES[ rnode ] );
      p_EQUATION -> y() = x_midpt;
      // Update the equation to the mid point position
      p_EQUATION -> update( temp_dummy );
      // loop over all the variables and fill the matrix
      for ( unsigned var = 0; var < order; ++var )
      {
        std::size_t placement_row( row + var );
        // deriv at the MID POINT between nodes
        ( *p_A ) ( placement_row, order * rnode + var ) = invh;
        ( *p_A ) ( placement_row, order * lnode + var ) = -invh;
        // add the Jacobian terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )
        {
          ( *p_A ) ( placement_row, order * lnode + i ) -= 0.5 * p_EQUATION -> jacobian()( var, i );
          ( *p_A ) ( placement_row, order * rnode + i ) -= 0.5 * p_EQUATION -> jacobian()( var, i );
        }
        // RHS
        for ( unsigned i = 0; i < order; ++i )
        {
          ( *p_B ) ( placement_row, order * lnode + i ) -= 0.5 * p_EQUATION -> mass()( var, i );
          ( *p_B ) ( placement_row, order * rnode + i ) -= 0.5 * p_EQUATION -> mass()( var, i );
        }
      }
      // increment the row
      row += order;
    }

    // update the BC residuals for the current iteration
    p_RIGHT_RESIDUAL -> update( temp_dummy );
    // add the (linearised) LHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_RIGHT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        ( *p_A ) ( row, order * ( num_of_nodes - 1 ) + var ) = p_RIGHT_RESIDUAL -> jacobian()( i, var );
        ( *p_B ) ( row, order * ( num_of_nodes - 1 ) + var ) = 0.0;
      }
      ++row;
    }

    if ( row != num_of_nodes * order )
    {
      std::string problem( "\n The ODE_BVP has an incorrect number of boundary conditions. \n" );
      throw ExceptionRuntime( problem );
    }

    // Dodgy hack alert: rather than flipping the signs in the above formulation, we will flip
    // them here to ensure that the RHS matrix (probably?) has the semi-definite form
    // required by the ARPACK solver.
    ( *p_A ).scale( -1.0 );
    ( *p_B ).scale( -1.0 );

  }


  // the templated versions that we require are:
  template class ODE_EVP<double>
  ;
  template class ODE_EVP<std::complex<double> >
  ;


} // end namespace

