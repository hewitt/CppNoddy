/// \file PDE_IBVP.cpp
/// A specification of a class for an \f$ n^{th} \f$-order IBVP of the form
/// \f[ M( {\underline f}(y,t), y, t )\cdot {\underline f}_t (y,t)+ {\underline f}_y (y,t) = {\underline R}( {\underline f}(y,t), y, t )\,, \f]
/// subject to \f$ n \f$ conditions defined at \f$ y = y_{left} \f$ and
/// \f$ y_{right} \f$ for some components of \f$ {\underline f}(y) \f$.
/// Here \f$ M \f$ is a mass matrix.
/// The solution at the new time step \f$ t+\Delta t \f$ is
/// \f[ {\underline f}^{new} = {\underline F} + {\underline g} \f]
/// where \f$ {\underline F} \f$ is the current guess at the solution
/// and \f$ {\underline g} \f$ is the linearised correction.
/// The solution at the previous time \f$ t \f$ is
/// \f[ {\underline f}^{old} = {\underline O} \f]
/// A Crank-Nicolson method is employed with the linearised problem at the mid-time point
/// \f$ t + \Delta t /2 \f$ being:
/// \f[ \frac{2}{\Delta t} M \cdot {\underline g } + {\underline g}_y - J \cdot {\underline g} + J^{(M)} \cdot \frac{\underline F - \underline O}{\Delta t} \cdot {\underline g}  = 2 {\underline R} - ({\underline F}_y + {\underline O}_y ) - \frac{2}{\Delta t} M \cdot ( {\underline F} - {\underline O} ) \f]
/// Where \f$ M, J, J^{(M)}, R \f$ are evaluated at the mid-time step with arguments \f$ \left ( \frac{\underline F + \underline O}{2}, y, t + \frac{\Delta t}{2} \right ) \f$,
/// with \f$ J^{(M)} \f$ denoting the Jacobian of the mass matrix \f$ \partial M_{ij} / \partial f_k \f$.
/// This problem is solved by second-order central differencing the equation at
/// the spatial (\f$ y \f$) inter-node mid points.

#include <string>
#include <cassert>

#include <PDE_IBVP.h>
#include <Equation_with_mass.h>
#include <Residual_with_coords.h>
#include <Exceptions.h>
#include <BandedLinearSystem.h>
#include <Timer.h>
#include <Utility.h>

namespace CppNoddy
{

  template <typename _Type>
  PDE_IBVP<_Type>::PDE_IBVP( Equation_with_mass<_Type >* ptr_to_equation,
                             const DenseVector<double> &nodes,
                             Residual_with_coords<_Type>* ptr_to_left_residual,
                             Residual_with_coords<_Type>* ptr_to_right_residual ) :
      TOL( 1.e-8 ),
      T( 0.0 ),
      MAX_ITERATIONS( 12 ),
      p_EQUATION( ptr_to_equation ),
      p_LEFT_RESIDUAL( ptr_to_left_residual ),
      p_RIGHT_RESIDUAL( ptr_to_right_residual )
  {
    SOLN = OneD_Node_Mesh<_Type>( nodes, p_EQUATION -> get_order() );
    PREV_SOLN = SOLN;
    p_EQUATION -> t() = T;
#ifdef TIME
    // timers
    T_ASSEMBLE = Timer( "Assembling of the matrix (incl. equation updates):" );
    T_SOLVE = Timer( "Solving of the matrix:" );
#endif
  }

  template <typename _Type>
  PDE_IBVP<_Type>::~PDE_IBVP()
  {
#ifdef TIME
    std::cout << "\n";
    T_ASSEMBLE.stop();
    T_ASSEMBLE.print();
    T_SOLVE.stop();
    T_SOLVE.print();
#endif
  }

  template <typename _Type>
  void PDE_IBVP<_Type>::step2( const double& dt )
  {
    // the order of the problem
    unsigned order( p_EQUATION -> get_order() );
    // get the number of nodes in the mesh
    // -- this may have been refined by the user since the last call.
    unsigned ny( SOLN.get_nnodes() );
    // measure of maximum residual
    double max_residual( 1.0 );
    // iteration counter
    int counter = 0;
    // store this soln as the 'previous SOLUTION'
    PREV_SOLN = SOLN;
    // ANY LARGE STORAGE USED IN THE MAIN LOOP IS
    // DEFINED HERE TO AVOID REPEATED CONSTRUCTION.
    // Note we blank the A matrix after every iteration.
    //
    // Banded LHS matrix - max obove diagonal band width is
    // from first variable at node i to last variable at node i+1
    BandedMatrix<_Type> a( ny * order, 2 * order - 1, 0.0 );
    // RHS
    DenseVector<_Type> b( ny * order, 0.0 );
    // linear solver definition
#ifdef LAPACK
    BandedLinearSystem<_Type> system( &a, &b, "lapack" );
#else
    BandedLinearSystem<_Type> system( &a, &b, "native" );
#endif
    // loop until converged or too many iterations
    do
    {
      // iteration counter
      ++counter;
#ifdef TIME
      T_ASSEMBLE.start();
#endif
      assemble_matrix_problem( a, b, dt );
      max_residual = b.inf_norm();
#ifdef DEBUG
      std::cout << " PDE_IBVP.solve : Residual_max = " << max_residual << " tol = " << TOL << "\n";
#endif
#ifdef TIME
      T_ASSEMBLE.stop();
      T_SOLVE.start();
#endif
      // linear solver
      system.solve();
      // keep the solution in a OneD_GenMesh object
      for ( std::size_t var = 0; var < order; ++var )
      {
        for ( std::size_t i = 0; i < ny; ++i )
        {
          SOLN( i, var ) += b[ i * order + var ];
        }
      }
#ifdef TIME
      T_SOLVE.stop();
#endif
    }
    while ( ( max_residual > TOL ) && ( counter < MAX_ITERATIONS ) );
    if ( counter >= MAX_ITERATIONS )
    {
      std::string problem( "\n The PDE_IBVP.step2 method took too many iterations. \n" );
      throw ExceptionItn( problem, counter, max_residual );
    }
#ifdef DEBUG
    std::cout << "[DEBUG] time-like variable = " << T << "\n";
#endif
    // set the time to the updated level
    T += dt;
  }


  template <typename _Type>
  void PDE_IBVP<_Type>::assemble_matrix_problem( BandedMatrix<_Type>& a, DenseVector<_Type>& b, const double& dt )
  {
    // clear the Jacobian matrix
    a.assign( 0.0 );
    // inverse of the time step
    const double inv_dt( 1. / dt );
    // the order of the problem
    const unsigned order( p_EQUATION -> get_order() );
    // number of spatial nodes
    const unsigned ny( SOLN.get_nnodes() );
    // row counter
    std::size_t row( 0 );
    // a matrix that is used in the Jacobian of the mass matrix terms
    DenseMatrix<_Type> h( order, order, 0.0 );
    // local state variable and functions
    DenseVector<_Type> F_midpt( order, 0.0 );
    DenseVector<_Type> O_midpt( order, 0.0 );
    DenseVector<_Type> state( order, 0.0 );
    DenseVector<_Type> state_dt( order, 0.0 );
    DenseVector<_Type> state_dy( order, 0.0 );
    // BCn equation is evaluated at the next time step
    p_LEFT_RESIDUAL -> coord( 0 ) = T + dt;
    // update the BC residuals for the current iteration
    p_LEFT_RESIDUAL -> update( SOLN.get_nodes_vars( 0 ) );
    // add the (linearised) LHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_LEFT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, var ) = p_LEFT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - p_LEFT_RESIDUAL -> residual()[ i ];
      ++row;
    }
    // inner nodes of the mesh, node = 0,1,2,...,N-2
    for ( std::size_t node = 0; node <= ny - 2; ++node )
    {
      const std::size_t l_node = node;
      const std::size_t r_node = node + 1;
      // inverse of step size
      const double inv_dy = 1. / ( SOLN.coord( r_node ) - SOLN.coord( l_node ) );
      // set the current solution at this node by 2nd order evaluation at mid point
      for ( unsigned var = 0; var < order; ++var )
      {
        const _Type F_midpt = ( SOLN( l_node, var ) + SOLN( r_node, var ) ) / 2.;
        const _Type O_midpt = ( PREV_SOLN( l_node, var ) + PREV_SOLN( r_node, var ) ) / 2.;
        state_dy[ var ] = ( SOLN( r_node, var ) - SOLN( l_node, var )
                            + PREV_SOLN( r_node, var ) - PREV_SOLN( l_node, var ) ) * inv_dy / 2.;
        state[ var ] = ( F_midpt + O_midpt ) / 2.;
        state_dt[ var ] = ( F_midpt - O_midpt ) * inv_dt;
      }
      // set the equation's y & t values to be mid points
      p_EQUATION -> y() = 0.5 * ( SOLN.coord( l_node ) + SOLN.coord( r_node ) );
      p_EQUATION -> t() = T + dt / 2;
      // Update the equation to the mid point position
      p_EQUATION -> update( state );
      // evaluate the Jacobian of mass contribution multiplied by state_dt
      p_EQUATION -> get_jacobian_of_mass_mult_vector( state, state_dt, h );
      // loop over all the variables
      //
      // to avoid repeated mapping arithmetic within operator() of the
      // BandedMatrix class we'll access the matrix with the iterator.
      //
      // il = position for (0, lnode*order)
      typename BandedMatrix<_Type>::elt_iter l_iter( a.get_elt_iter( 0, l_node * order ) );
      // ir = position for (0, rnode*order)
      typename BandedMatrix<_Type>::elt_iter r_iter( a.get_elt_iter( 0, r_node * order ) );
      for ( unsigned var = 0; var < order; ++var )
      {
        // offset for (r,c) -> (r+row, c+var)
        const std::size_t offset( var * a.noffdiag() * 3 + row );
        // deriv at the MID POINT between nodes
        //a( row, order * l_node + var ) = -inv_dy;
        *( l_iter + offset ) = -inv_dy;
        //a( row, order * r_node + var ) = inv_dy;
        *( r_iter + offset ) = inv_dy;
        // add the matrix mult terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )  // dummy index
        {
          // add the Jacobian terms
          // indirect access: a( row, order * l_node + i ) -= jac_midpt( var, i ) * 0.5;
          const std::size_t offset( i * a.noffdiag() * 3 + row );
          typename BandedMatrix<_Type>::elt_iter left( l_iter + offset );
          typename BandedMatrix<_Type>::elt_iter right( r_iter + offset );
          *left -= p_EQUATION -> jacobian()( var, i ) * 0.5;
          // indirect access: a( row, order * r_node + i ) -= jac_midpt( var, i ) * 0.5;
          *right -= p_EQUATION -> jacobian()( var, i ) * 0.5;
          // add the Jacobian of mass terms
          // indirect access: a( row, order * l_node + i ) += h( var, i ) * 0.5;
          *left += h( var, i ) * 0.5;
          // indirect access: a( row, order * r_node + i ) += h( var, i ) * 0.5;
          *right += h( var, i ) * 0.5;
          // add the mass matrix terms
          // indirect access: a( row, order * l_node + i ) += mass_midpt( var, i ) * inv_dt;
          *left += p_EQUATION -> mass()( var, i ) * inv_dt;
          // indirect access: a( row, order * r_node + i ) += mass_midpt( var, i ) * inv_dt;
          *right += p_EQUATION -> mass()( var, i ) * inv_dt;
        }
        // RHS
        b[ row ] = p_EQUATION -> residual()[ var ] - state_dy[ var ];
        b[ row ] -= Utility::dot( p_EQUATION -> mass()[ var ], state_dt );
        b[ row ] *= 2;
        // increment the row
        row += 1;
      }
    }
    // BCn equation is evaluated at the next step time point as
    // they cannot depend on d/dt terms
    p_RIGHT_RESIDUAL -> coord( 0 ) = T + dt;
    // update the BC residuals for the current iteration
    p_RIGHT_RESIDUAL -> update( SOLN.get_nodes_vars( ny - 1 ) );
    // add the (linearised) RHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_RIGHT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at RHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, order * ( ny - 1 ) + var ) = p_RIGHT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - p_RIGHT_RESIDUAL -> residual()[ i ];
      ++row;
    }
#ifdef PARANOID
    if ( row != ny * order )
    {
      std::string problem( "\n The ODE_BVP has an incorrect number of boundary conditions. \n" );
      throw ExceptionRuntime( problem );
    }
#endif
  }

  // the templated versions that we require are:
  template class PDE_IBVP<double>
  ;
  template class PDE_IBVP< D_complex >
  ;

} // end namespace

