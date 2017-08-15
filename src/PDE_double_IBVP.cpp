/// \file PDE_double_IBVP.cpp


#include <string>
#include <cassert>

#include <PDE_double_IBVP.h>
#include <BandedLinearSystem.h>
#include <Equation_3matrix.h>
#include <Exceptions.h>
#include <Utility.h>
#include <Timer.h>

//#define DEBUG

namespace CppNoddy
{
  template <typename _Type>
  PDE_double_IBVP<_Type>::PDE_double_IBVP(
    Equation_3matrix<_Type >* ptr_to_equation,
    const DenseVector<double> &xnodes,
    const DenseVector<double> &ynodes,
    Residual_with_coords<_Type>* ptr_to_bottom_residual,
    Residual_with_coords<_Type>* ptr_to_top_residual ) :
      TOL( 1.e-8 ),
      T( 0.0 ),
      MAX_ITERATIONS( 20 ),
      p_EQUATION( ptr_to_equation ),
      p_BOTTOM_RESIDUAL( ptr_to_bottom_residual ),
      p_TOP_RESIDUAL( ptr_to_top_residual ),
      UPDATED( false )
  {
    // create the 2D mesh for the current time level and previous time level
    SOLN = TwoD_Node_Mesh<_Type>( xnodes, ynodes, p_EQUATION -> get_order() );
    // copy construct the previous solution's storage
    PREV_SOLN = SOLN;
    // initialise the eqn object
    // "x"
    p_EQUATION -> coord(2) = xnodes[0];
#ifdef TIME
    // timers
    T_ASSEMBLE = Timer( "Assembling of the matrix (incl. equation updates):" );
    T_SOLVE = Timer( "Solving of the matrix:" );
#endif
  }

  template <typename _Type>
  PDE_double_IBVP<_Type>::~PDE_double_IBVP()
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
  void PDE_double_IBVP<_Type>::step2( const double& dt )
  {
    DT = dt;
    if ( !UPDATED )
    {
      std::string problem;
      problem = " You have called the PDE_double_IBVP::step2 method without calling\n";
      problem += " the PDE_double_IBVP::update_previous_solution method first. This\n";
      problem += " method is required to copy/store the previous time level's solution\n";
      problem += " and then allow you to specify/alter the x-initial condition by\n";
      problem += " for the current time level.";
      throw ExceptionRuntime( problem );
    }
    // the order of the problem
    unsigned order( p_EQUATION -> get_order() );
    // get the number of nodes in the mesh
    // -- this may have been refined by the user since the last call.
    unsigned nx( SOLN.get_nnodes().first );
    unsigned ny( SOLN.get_nnodes().second );
    // measure of maximum residual
    double max_residual( 1.0 );
    // the current solution moves to the previous solution
    // ANY LARGE STORAGE USED IN THE MAIN LOOP IS
    // DEFINED HERE TO AVOID REPEATED CONSTRUCTION.
    // Note we blank the A matrix after every iteration.
    //
    // Banded LHS matrix - max obove diagonal band wTidth is
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
    for ( std::size_t j = 0; j < nx - 1; ++j )
    {
      //// next x-soln guess is the previous x-soln
      //for ( std::size_t var = 0; var < order; ++var )
      //{
        //for ( std::size_t i = 0; i < ny; ++i )
        //{
          //SOLN( j + 1, i, var ) = SOLN( j, i, var );
        //}
      //}
      //
      // iteration counter
      int counter = 0;
      // loop until converged or too many iterations
      do
      {
        // iteration counter
        ++counter;
        // time the assemble phase
#ifdef TIME
        T_ASSEMBLE.start();
#endif
        assemble_matrix_problem( a, b, j );
        max_residual = b.inf_norm();
#ifdef DEBUG
        std::cout.precision( 12 );
        std::cout << " PDE_double_IBVP.step2 : x_j = " << SOLN.coord(j,0).first << " Residual_max = " << max_residual << " tol = " << TOL << " counter = " << counter << "\n";
#endif
#ifdef TIME
        T_ASSEMBLE.stop();
        T_SOLVE.start();
#endif
        // linear solver
        system.solve();
        // keep the solution in a OneD_GenMesh object
        for ( std::size_t i = 0; i < ny; ++i )
        {
          for ( std::size_t var = 0; var < order; ++var )
          {
            SOLN( j + 1, i, var ) += b[ i * order + var ];
          }
        }
#ifdef TIME
        T_SOLVE.stop();
#endif
      }
      while ( ( max_residual > TOL ) && ( counter < MAX_ITERATIONS ) );
      if ( counter >= MAX_ITERATIONS )
      {
        std::string problem( "\n The PDE_double_IBVP.step2 method took too many iterations. \n" );
        throw ExceptionItn( problem, counter, max_residual );
      }
    }
    // set the time to the updated level
    T += DT;
#ifdef DEBUG
    std::cout << "[DEBUG] solved at t = " << T << "\n";
#endif
    UPDATED = false;
  }


  template <typename _Type>
  void PDE_double_IBVP<_Type>::assemble_matrix_problem(
    BandedMatrix<_Type>& a,
    DenseVector<_Type>& b,
    const std::size_t& j )
  {
    // clear the Jacobian matrix
    a.assign( 0.0 );
    // inverse of the t-step
    const double inv_dt( 1. / DT );
    // inverse of the x-step
    const double inv_dx( 1. / ( SOLN.coord( j + 1, 0 ).first - SOLN.coord( j, 0 ).first ) );
    // mid point of the x variable
    const double x_midpt = 0.5 * ( SOLN.coord( j + 1, 0 ).first + SOLN.coord( j, 0 ).first );
    // the order of the problem
    const unsigned order( p_EQUATION -> get_order() );
    // number of spatial nodes in the y-direction
    const unsigned ny( SOLN.get_nnodes().second );
    // row counter
    std::size_t row( 0 );
    // local state variables
    DenseVector<_Type> state( order, 0.0 );
    DenseVector<_Type> state_dx( order, 0.0 );
    DenseVector<_Type> state_dt( order, 0.0 );
    DenseVector<_Type> state_dy( order, 0.0 );
    // Current and old state
    DenseVector<_Type> F_midpt( order, 0.0 );
    DenseVector<_Type> O_midpt( order, 0.0 );
    // these store the Jacobian of the mass matrices times a state vector
    DenseMatrix<_Type> h0( order, order, 0.0 );
    DenseMatrix<_Type> h1( order, order, 0.0 );
    DenseMatrix<_Type> h2( order, order, 0.0 );
    // set the x-t coordinates in the BC
    p_BOTTOM_RESIDUAL -> coord( 0 ) = SOLN.coord( j + 1, 0 ).first;
    p_BOTTOM_RESIDUAL -> coord( 1 ) = T + DT;
    // update the BC residuals for the current iteration
    p_BOTTOM_RESIDUAL -> update( SOLN.get_nodes_vars( j + 1, 0 ) );
    // add the (linearised) bottom BCs to the matrix problem
    for ( unsigned i = 0; i < p_BOTTOM_RESIDUAL-> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, var ) = p_BOTTOM_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = -p_BOTTOM_RESIDUAL -> residual()[ i ];
      ++row;
    }
    // inner nodes of the mesh, node = 0,1,2,...,ny-2
    for ( std::size_t node = 0; node <= ny - 2; ++node )
    {
      // bottom and top nodes of the mid-y-point we're considering
      const std::size_t b_node = node;
      const std::size_t t_node = node + 1;
      //
      const double inv_dy = 1. / ( SOLN.coord( j, t_node ).second - SOLN.coord( j, b_node ).second );
      // work out state, state_dx, state_dt, state_dy ... all at the mid-point in x,y,t.
      for ( unsigned var = 0; var < order; ++var )
      {
        const double Fj_midpt = ( SOLN( j, b_node, var ) + SOLN( j, t_node, var ) ) / 2;
        const double Oj_midpt = ( PREV_SOLN( j, b_node, var ) + PREV_SOLN( j, t_node, var ) ) / 2;
        const double Fjp1_midpt = ( SOLN( j + 1, b_node, var ) + SOLN( j + 1, t_node, var ) ) / 2;
        const double Ojp1_midpt = ( PREV_SOLN( j + 1, b_node, var ) + PREV_SOLN( j + 1, t_node, var ) ) / 2;
        state[ var ] = ( Fj_midpt + Fjp1_midpt + Oj_midpt + Ojp1_midpt ) / 4;
        state_dx[ var ] = ( Fjp1_midpt - Fj_midpt + Ojp1_midpt - Oj_midpt ) * inv_dx / 2;
        state_dt[ var ] = ( Fjp1_midpt + Fj_midpt - Ojp1_midpt - Oj_midpt ) * inv_dt / 2;
        state_dy[ var ] = ( SOLN( j + 1, t_node, var ) - SOLN( j + 1, b_node, var )
                            + SOLN( j, t_node, var ) - SOLN( j, b_node, var )
                            + PREV_SOLN( j + 1, t_node, var ) - PREV_SOLN( j + 1, b_node, var )
                            + PREV_SOLN( j, t_node, var ) - PREV_SOLN( j, b_node, var ) ) * inv_dy / 4;
      }
      //
      double y_midpt = 0.5 * ( SOLN.coord( j, b_node ).second + SOLN.coord( j, t_node ).second );
      // set the coord in the equation object
      p_EQUATION -> coord(0) = y_midpt;
      p_EQUATION -> coord(1) = T + DT / 2.;
      p_EQUATION -> coord(2) = x_midpt;
      // Update the equation to the mid point position
      p_EQUATION -> update( state );
      // evaluate the Jacobian of mass contribution multiplied by state_dy
      p_EQUATION -> get_jacobian_of_matrix0_mult_vector( state, state_dy, h0 );
      // evaluate the Jacobian of mass contribution multiplied by state_dt
      p_EQUATION -> get_jacobian_of_matrix1_mult_vector( state, state_dt, h1 );
      // evaluate the Jacobian of mass contribution multiplied by state_dx
      p_EQUATION -> get_jacobian_of_matrix2_mult_vector( state, state_dx, h2 );
      // loop over all the variables
      // il = position for (0, b_node*order)
      typename BandedMatrix<_Type>::elt_iter b_iter( a.get_elt_iter( 0, b_node * order ) );
      // ir = position for (0, t_node*order)
      typename BandedMatrix<_Type>::elt_iter t_iter( a.get_elt_iter( 0, t_node * order ) );
      for ( unsigned var = 0; var < order; ++var )
      {
        // add the matrix mult terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )  // dummy index
        {
          const std::size_t offset( i * a.noffdiag() * 3 + row );
          typename BandedMatrix<_Type>::elt_iter bottom( b_iter + offset );
          typename BandedMatrix<_Type>::elt_iter top( t_iter + offset );
          // add the Jacobian terms
          //a( row, order * b_node + i ) -=  p_EQUATION -> jacobian()( var, i ) / 2;
          *bottom -= p_EQUATION -> jacobian()( var, i ) / 2;
          //a( row, order * t_node + i ) -=  p_EQUATION -> jacobian()( var, i ) / 2;
          *top -= p_EQUATION -> jacobian()( var, i ) / 2;
          // add the Jacobian of matrix terms
          // indirect access: a( row, order * l_node + i ) += h1( var, i ) * 0.5;
          *bottom += h0( var, i ) * 0.5;
          // indirect access: a( row, order * r_node + i ) += h1( var, i ) * 0.5;
          *top += h0( var, i ) * 0.5;
          //a( row, order * b_node + i ) += h1( var, i ) * .5;
          *bottom += h1( var, i ) * .5;
          //a( row, order * t_node + i ) += h1( var, i ) * .5;
          *top += h1( var, i ) * .5;
          // add the Jacobian of matrix terms
          //a( row, order * b_node + i ) += h2( var, i ) * .5;
          *bottom += h2( var, i ) * .5;
          //a( row, order * t_node + i ) += h2( var, i ) * .5;
          *top += h2( var, i ) * .5;
          // add the matrix terms
          
          // indirect access: a( row, order * l_node + i ) -= mass_midpt( var, i ) * inv_dy;
          *bottom -= p_EQUATION -> matrix0()( var, i ) * inv_dy;
          // indirect access: a( row, order * r_node + i ) += mass_midpt( var, i ) * inv_dy;
          *top += p_EQUATION -> matrix0()( var, i ) * inv_dy;
          
          //a( row, order * b_node + i ) += p_EQUATION -> mass2()( var, i ) * inv_dt;
          *bottom += p_EQUATION -> matrix1()( var, i ) * inv_dt;
          //a( row, order * t_node + i ) += p_EQUATION -> mass2()( var, i ) * inv_dt;
          *top += p_EQUATION -> matrix1()( var, i ) * inv_dt;
          //a( row, order * b_node + i ) += p_EQUATION -> mass1()( var, i ) * inv_dx;
          *bottom += p_EQUATION -> matrix2()( var, i ) * inv_dx;
          //a( row, order * t_node + i ) += p_EQUATION -> mass1()( var, i ) * inv_dx;
          *top += p_EQUATION -> matrix2()( var, i ) * inv_dx;
        }
        // RHS
        b[ row ] = p_EQUATION -> residual()[ var ];
        b[ row ] -= Utility::dot( p_EQUATION -> matrix0()[ var ], state_dy );
        b[ row ] -= Utility::dot( p_EQUATION -> matrix1()[ var ], state_dt );
        b[ row ] -= Utility::dot( p_EQUATION -> matrix2()[ var ], state_dx );
        b[ row ] *= 4;
        // increment the row
        row += 1;
      }
    }
    // set the x-t coordinates in the BC
    p_TOP_RESIDUAL -> coord( 0 ) = SOLN.coord( j + 1, ny - 1 ).first;
    p_TOP_RESIDUAL -> coord( 1 ) = T + DT;
    // update the BC residuals for the current iteration
    p_TOP_RESIDUAL -> update( SOLN.get_nodes_vars( j + 1, ny - 1 ) );
    // add the (linearised) RHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_TOP_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at RHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, order * ( ny - 1 ) + var ) = p_TOP_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - p_TOP_RESIDUAL -> residual()[ i ];
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
  template class PDE_double_IBVP<double>
  ;

} // end namespace

