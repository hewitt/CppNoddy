/// \file ODE_BVP.cpp
/// An implementation of a class for an \f$ n^{th} \f$-order ODE BVP defined by
/// \f[ {\underline f}^\prime (x) = {\underline R}( {\underline f}(x), x )\,, \f]
/// subject to \f$ n \f$ Dirichlet conditions defined at \f$ x = x_{left} \f$ or
/// \f$ x_{right} \f$ for some components of \f$ {\underline f}(x) \f$. The system
/// is solved by applying Newton iteration, with the intermediate problem:
/// \f[ {\underline g}^\prime (x) - \frac{\partial {\underline R}}{\partial \underline f} \Big \vert_{\underline F} \,\,{\underline g}(x) = {\underline R}( {\underline F}(x), x) - {\underline F}^\prime (x) \,, \f]
/// for the corrections \f$ \underline g(x) \f$ to the current approximation
/// to the solution \f$ {\underline F}(x) \f$.The numerical scheme can be
/// run for any given distribution of nodes and can adapt the nodal
/// positions based on residual evaluations (including both refinement and unrefinement).

#include <string>
#include <cassert>
#include <utility>
#include <iostream>
#include <typeinfo>

#include <ODE_BVP.h>
#include <Residual_with_coords.h>
#include <ArcLength_base.h>
#include <BandedLinearSystem.h>
#include <Residual.h>
#include <Utility.h>
#include <Exceptions.h>
#include <Timer.h>
#include <OneD_Node_Mesh.h>
#include <Equation_with_mass.h>

namespace CppNoddy
{

  template <typename _Type, typename _Xtype>
  ODE_BVP<_Type, _Xtype>::ODE_BVP( Residual_with_coords<_Type, _Xtype>* ptr_to_equation,
                                   const DenseVector<_Xtype> &nodes,
                                   Residual<_Type>* ptr_to_left_residual,
                                   Residual<_Type>* ptr_to_right_residual ) :
      ArcLength_base<_Type> (),
      MAX_ITERATIONS( 12 ),
      TOL( 1.e-8 ),
      LAST_DET_SIGN( 0 ),
      p_EQUATION( ptr_to_equation ),
      p_LEFT_RESIDUAL( ptr_to_left_residual ),
      p_RIGHT_RESIDUAL( ptr_to_right_residual ),
      MONITOR_DET( true )
  {
    // set up the solution mesh using these nodes
    SOLUTION = OneD_Node_Mesh<_Type, _Xtype>( nodes, p_EQUATION -> get_order() );
    if ( ( p_LEFT_RESIDUAL -> get_number_of_vars() != p_EQUATION -> get_order() ) ||
         ( p_RIGHT_RESIDUAL -> get_number_of_vars() != p_EQUATION -> get_order() ) ||
         ( p_LEFT_RESIDUAL -> get_order() + p_RIGHT_RESIDUAL -> get_order() != p_EQUATION -> get_order() ) )
    {
      std::string problem;
      problem = "It looks like the ODE_BVP equation and boundary conditions are\n";
      problem += "not well posed. The number of variables for each boundary condition\n";
      problem += "has to be the same as the order of the equation. Also the order of\n";
      problem += "both boundary conditions has to sum to the order of the equation.\n";
      throw ExceptionBifurcation( problem );
    }

#ifdef TIME
    // timers
    T_ASSEMBLE = Timer( "ODE_BVP: Assembling of the matrix (incl. equation updates):" );
    T_SOLVE = Timer( "ODE_BVP: Solving of the matrix:" );
    T_REFINE = Timer( "ODE_BVP: Refining the mesh:" );
#endif
  }

  template <typename _Type, typename _Xtype>
  ODE_BVP<_Type, _Xtype>::~ODE_BVP()
  {
#ifdef TIME
    std::cout << "\n";
    T_ASSEMBLE.stop();
    T_ASSEMBLE.print();
    T_SOLVE.stop();
    T_SOLVE.print();
    T_REFINE.stop();
    T_REFINE.print();
#endif
  }


  template <typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::solve2()
  {
    // the order of the problem
    unsigned order( p_EQUATION -> get_order() );
    // get the number of nodes in the mesh
    // -- this may have been refined by the user since the
    // last call.
    unsigned n( SOLUTION.get_nnodes() );
    // measure of maximum residual
    double max_residual( 1.0 );
    // iteration counter
    int counter = 0;
    // determinant monitor
    int det_sign( LAST_DET_SIGN );

    // ANY LARGE STORAGE USED IN THE MAIN LOOP IS
    // DEFINED HERE TO AVOID REPEATED CONSTRUCTION.
    // Note we blank the A matrix after every iteration.
    //
    // Banded LHS matrix - max obove diagonal band width is
    // from first variable at node i to last variable at node i+1
    BandedMatrix<_Type> a( n * order, 2 * order - 1, 0.0 );
    // RHS
    DenseVector<_Type> b( n * order, 0.0 );
    // linear solver definition
#ifdef LAPACK
    BandedLinearSystem<_Type> system( &a, &b, "lapack" );
#else
    BandedLinearSystem<_Type> system( &a, &b, "native" );
#endif
    // pass the local monitor_det flag to the linear system
    system.set_monitor_det( MONITOR_DET );

    // loop until converged or too many iterations
    do
    {
      // iteration counter
      ++counter;
#ifdef TIME
      T_ASSEMBLE.stop();
#endif
      assemble_matrix_problem( a, b );
      actions_before_linear_solve( a, b );
      max_residual = b.inf_norm();
#ifdef DEBUG
      std::cout << " ODE_BVP.solve : Residual_max = " << max_residual << " tol = " << TOL << "\n";
#endif
#ifdef TIME
      T_ASSEMBLE.start();
      T_SOLVE.start();
#endif
      // linear solver
      system.solve();
      // keep the solution in a OneD_Node_Mesh object
      for ( std::size_t var = 0; var < order; ++var )
      {
        for ( std::size_t i = 0; i < n; ++i )
        {
          SOLUTION( i, var ) += b[ i * order + var ];
        }
      }
#ifdef TIME
      T_SOLVE.stop();
#endif
    }
    while ( ( max_residual > TOL ) && ( counter < MAX_ITERATIONS ) );
    if ( counter >= MAX_ITERATIONS )
    {
      std::string problem( "\n The ODE_BVP.solve2 method took too many iterations. \n" );
      throw ExceptionItn( problem, counter, max_residual );
    }

    // last_det_sign is set to be 0 in ctor, it must be +/-1 when calculated
    if ( MONITOR_DET )
    {
      det_sign = system.get_det_sign();
      if ( det_sign * LAST_DET_SIGN < 0 )
      {
        LAST_DET_SIGN = det_sign;
        std::string problem;
        problem = "[ INFO ] : Determinant monitor has changed signs in ODE_BVP.\n";
        problem += "[ INFO ] : Bifurcation detected.\n";
        throw ExceptionBifurcation( problem );
      }
      LAST_DET_SIGN = det_sign;
    }
  }

  template <typename _Type, typename _Xtype>
  double ODE_BVP<_Type, _Xtype>::arclength_solve( const double& step )
  {
    this -> ds() = step;
    // order of the equation
    unsigned order( p_EQUATION -> get_order() );
    // the number of nodes in the BVP
    unsigned n( SOLUTION.get_nnodes() );
    // Banded Jacobian
    BandedMatrix<_Type> Jac( n * order, 2 * order - 1, 0.0 );
    // residuals over all nodes vectors
    DenseVector<_Type> Res1( n * order, 0.0 );
    DenseVector<_Type> Res2( n * order, 0.0 );
    // RHS vectors for the linear solver(s)
    DenseVector<_Type> y( n * order, 0.0 );
    DenseVector<_Type> z( n * order, 0.0 );
#ifdef LAPACK
    BandedLinearSystem<_Type> system1( &Jac, &y, "lapack" );
    BandedLinearSystem<_Type> system2( &Jac, &z, "lapack" );
#else
    BandedLinearSystem<_Type> system1( &Jac, &y, "native" );
    BandedLinearSystem<_Type> system2( &Jac, &z, "native" );
#endif
    // we use the member data 'monitor_det' to set the linear system determinant monitor
    system1.set_monitor_det( MONITOR_DET );
    // make backups in case we can't find a converged solution
    DenseVector<_Type> backup_state( SOLUTION.vars_as_vector() );
    _Type backup_parameter( *( this -> p_PARAM ) );
    // we can generate a 1st-order guess for the next state and parameter
    DenseVector<_Type> x( this -> LAST_X + this -> X_DERIV_S * this -> DS );
    *( this -> p_PARAM ) = this -> LAST_PARAM + this -> PARAM_DERIV_S * this -> DS;
    // the class keeps the solution in a OneD_mesh object, so we update it here
    SOLUTION.set_vars_from_vector( x );
    // determinant monitor
    int det_sign( 0 );
    // a flag to indicate if we have made a successful step
    bool step_succeeded( false );
    // iteration counter
    int counter = 0;
    // loop until converged or too many iterations
    do
    {
      /// \todo y & z should be solved for simultaneously - but to do this
      /// we'd have to extend the LAPACK interface and/or provide the
      /// same feature for the native solver.
      ++counter;
      // extra constraint corresponding to the new unknow parameter (arclength)
      double E1 = this -> arclength_residual( x );
      // construct the Jacobian/residual matrix problem
      assemble_matrix_problem( Jac, Res1 );
      actions_before_linear_solve( Jac, Res1 );
      BandedMatrix<_Type> Jac_copy( Jac );
      y = Res1;
#ifdef DEBUG
      //std::cout << " [ DEBUG ] : max_residual = " << Res1.inf_norm() << "\n";
#endif
      if ( Res1.inf_norm() < TOL && counter > 1 )
      {
        step_succeeded = true;
        break;
      }
      try
      {
        system1.solve();
        det_sign = system1.get_det_sign();
      }
      catch ( ExceptionExternal )
      {
        step_succeeded = false;
        break;
      }
      // derivs w.r.t parameter
      const double delta( 1.e-8 );
      *( this -> p_PARAM ) += delta;
      // we actually just want the Res2 & e2 vectors, so we still
      // use the original Jacobian j ... but it was overwritten by the
      // previous solve.
      assemble_matrix_problem( Jac, Res2 );
      //Jac = Jac_copy;
      double E2 = this -> arclength_residual( x );
      actions_before_linear_solve( Jac, Res2 );
      *( this -> p_PARAM )  -= delta;
      DenseVector<_Type> dRes_dp(  ( Res2 - Res1 ) / delta );
      double dE_dp = ( E2 - E1 ) / delta;
      z = dRes_dp;
      try
      {
        system2.solve();
      }
      catch ( ExceptionExternal )
      {
        step_succeeded = false;
        break;
      }
      // this is the extra (full) row in the augmented matrix problem
      // bottom row formed from dE/dx_j
      DenseVector<_Type> Jac_E( this -> Jac_arclength_residual( x ) );
      // given the solutions y & z, the bordering algo provides the
      // solution to the full sparse system via
      _Type delta_p = - ( E1 + Utility::dot( Jac_E, y ) ) /
                      ( dE_dp + Utility::dot( Jac_E, z ) );
      DenseVector<_Type> delta_x = y + z * delta_p;
#ifdef DEBUG
      std::cout << " [ DEBUG ] : Arclength_solve, dp = " << delta_p
                << " dx.inf_norm() = " << delta_x.inf_norm()
                << " theta = " << this -> THETA << "\n";
#endif
      // update the state variables and the parameter with corrections
      x += delta_x;
      *( this -> p_PARAM ) += delta_p;
      // the corrections must be put back into the mesh container
      SOLUTION.set_vars_from_vector( x );
      // converged?
      if ( delta_x.inf_norm() < TOL )
      {
        step_succeeded = true;
        break;
      }
      // too many iterations?
      if ( counter > MAX_ITERATIONS )
      {
        step_succeeded = false;
        break;
      }
    }
    while ( true );
    //
    if ( !step_succeeded )
    {
      // if not a successful step then restore things
      SOLUTION.set_vars_from_vector( backup_state );
      *( this -> p_PARAM ) = backup_parameter;
      // reduce our step length
      this -> DS /= this -> ARCSTEP_MULTIPLIER;
#ifdef DEBUG
      std::cout << "[ DEBUG ] : REJECTING STEP \n";
      std::cout << "[ DEBUG ] : I decreased DS to " << this -> DS << "\n";
#endif
    }
    else
    {
      // update the variables needed for arc-length continuation
      update( SOLUTION.vars_as_vector() );
      if ( LAST_DET_SIGN * det_sign < 0 )
      {
        // a change in the sign of the determinant of the Jacobian
        // has been found
        LAST_DET_SIGN = det_sign;
        std::string problem;
        problem = "[ INFO ] : Determinant monitor has changed signs in the Newton class.\n";
        problem += "[ INFO ] : Bifurcation detected.\n";
        throw ExceptionBifurcation( problem );
      }

#ifdef DEBUG
      std::cout << " [ DEBUG ] : Number of iterations = " << counter << "\n";
      std::cout << " [ DEBUG ] : Parameter p = " << *( this -> p_PARAM )
                << "; arclength DS = " << this -> DS << "\n";
      std::cout << " [ DEBUG ] : Arclength THETA = " << this -> THETA << "\n";
#endif
      if ( counter > 8 || std::abs( this -> DS ) > this -> MAX_DS )
      {
        // converging too slowly, so decrease DS
        this -> DS /= this -> ARCSTEP_MULTIPLIER;
#ifdef DEBUG
        std::cout << " [ DEBUG ] : I decreased DS to " << this -> DS << "\n";
#endif
      }
      if ( counter < 4 )
      {
        if ( std::abs( this -> DS * this -> ARCSTEP_MULTIPLIER ) < this -> MAX_DS )
        {
          // converging too quickly, so increase DS
          this -> DS *= this -> ARCSTEP_MULTIPLIER;
#ifdef DEBUG
          std::cout << " [ DEBUG ] : I increased DS to " << this -> DS << "\n";
#endif
        }
      }
    }
    return this -> DS;
  }

  template <typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::solve( DenseVector<_Type>& state )
  {
    SOLUTION.set_vars_from_vector( state );
    try
    {
      solve2();
    }
    catch ( ExceptionBifurcation )
    {
      // dont throw from this routine, it's only needed for
      // the arclength_base class
    }
    state = SOLUTION.vars_as_vector();
  }


  template <typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::assemble_matrix_problem( BandedMatrix<_Type>& a, DenseVector<_Type>& b )
  {
    // clear the Jacobian matrix, since not all elts will be overwritten
    // in the routine below.
    a.assign( 0.0 );
    // the order of the problem
    unsigned order( p_EQUATION -> get_order() );
    // get the number of nodes in the mesh
    unsigned n( SOLUTION.get_nnodes() );
    // row counter
    std::size_t row( 0 );
    // Jacobian matrix for the equation
    DenseMatrix<_Type> Jac_midpt( order, order, 0.0 );
    // local state variable and functions
    DenseVector<_Type> F_midpt( order, 0.0 );
    DenseVector<_Type> R_midpt( order, 0.0 );
    // update the BC residuals for the current iteration
    p_LEFT_RESIDUAL -> update( SOLUTION.get_nodes_vars( 0 ) );
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
    for ( std::size_t node = 0; node <= n - 2; ++node )
    {
      const std::size_t lnode = node;
      const std::size_t rnode = node + 1;
      // inverse of step size
      const _Xtype invh = 1. / ( SOLUTION.coord( rnode ) - SOLUTION.coord( lnode ) );
      // set the current solution at this node by 2nd order evaluation at mid point
      for ( unsigned var = 0; var < order; ++var )
      {
        //F_midpt[ var ] = ( SOLUTION[ var ][ lnode ] + SOLUTION[ var ][ rnode ] ) / 2.;
        F_midpt[ var ] = ( SOLUTION( lnode, var ) + SOLUTION( rnode, var ) ) / 2.;
      }
      // mid point of the independent variable
      _Xtype y_midpt = ( SOLUTION.coord( lnode ) + SOLUTION.coord( rnode ) ) / 2.;
      // set the equation coord
      p_EQUATION -> y() = y_midpt;
      // Update the equation to the mid point position
      p_EQUATION -> update( F_midpt );
      // loop over all the variables
      for ( unsigned var = 0; var < order; ++var )
      {
        // deriv at the MID POINT between nodes
        a( row, order * rnode + var ) = invh;
        a( row, order * lnode + var ) = -invh;
        // add the Jacobian terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )
        {
          a( row, order * lnode + i ) -= p_EQUATION -> jacobian()( var, i ) / 2.;
          a( row, order * rnode + i ) -= p_EQUATION -> jacobian()( var, i ) / 2.;
        }
        // RHS
        b[ row ] = p_EQUATION -> residual()[ var ] - ( SOLUTION( rnode, var ) - SOLUTION( lnode, var ) ) * invh;
        // increment the row
        row += 1;
      }
    }
    // update the BC residuals for the current iteration
    p_RIGHT_RESIDUAL -> update( SOLUTION.get_nodes_vars( n - 1 ) );
    // add the (linearised) RHS BCs to the matrix problem
    for ( unsigned i = 0; i < p_RIGHT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at RHS of the domain
      for ( unsigned var = 0; var < order; ++var )
      {
        a( row, order * ( n - 1 ) + var ) = p_RIGHT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - p_RIGHT_RESIDUAL -> residual()[ i ];
      ++row;
    }
#ifdef PARANOID
    if ( row != n * order )
    {
      std::string problem( "\n The ODE_BVP has an incorrect number of boundary conditions. \n" );
      throw ExceptionRuntime( problem );
    }
#endif
  }


  // specialise to remove adaptivity from problems in the complex plane
  template<>
  std::pair< unsigned, unsigned > ODE_BVP<std::complex<double>, std::complex<double> >::adapt( const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt method has been called for a  \n";
    problem += " problem in the complex plane. This doesn't make sense \n";
    problem += " as the path is not geometrically defined. \n";
    throw ExceptionRuntime( problem );
  }

  // specialise to remove adaptivity from problems in the complex plane
  template<>
  std::pair< unsigned, unsigned > ODE_BVP<double, std::complex<double> >::adapt( const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt method has been called for a  \n";
    problem += " problem in the complex plane. This doesn't make sense \n";
    problem += " as the path is not geometrically defined. \n";
    throw ExceptionRuntime( problem );
  }

  // specialise to remove adaptivity from problems in the complex plane
  template<>
  void ODE_BVP<std::complex<double>, std::complex<double> >::adapt_until( const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt method has been called for a  \n";
    problem += " problem in the complex plane. This doesn't make sense \n";
    problem += " as the path is not geometrically defined. \n";
    throw ExceptionRuntime( problem );
  }

  // specialise to remove adaptivity from problems in the complex plane
  template<>
  void ODE_BVP<double, std::complex<double> >::adapt_until( const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt method has been called for a  \n";
    problem += " problem in the complex plane. This doesn't make sense \n";
    problem += " as the path is not geometrically defined. \n";
    throw ExceptionRuntime( problem );
  }

  template< typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::adapt_until( const double& adapt_tol )
  {
    std::pair< unsigned, unsigned > changes;
    do
    {
      changes = adapt( adapt_tol );
      solve2();
      std::cout << "[INFO] Adapting: " << changes.first << " " << changes.second << "\n";
    }
    while ( changes.first + changes.second != 0 );
  }

  template <typename _Type, typename _Xtype>
  std::pair< unsigned, unsigned > ODE_BVP<_Type, _Xtype>::adapt( const double& adapt_tol )
  {
#ifdef TIME
    T_REFINE.start();
#endif
    // the order of the problem
    unsigned order( p_EQUATION -> get_order() );
    // get the number of nodes in the mesh
    // -- this may have been refined by the user since the
    // last call.
    unsigned N( SOLUTION.get_nnodes() );
    // row counter
    std::size_t row( 0 );
    // local state variable and functions
    DenseVector<_Type> F_node( order, 0.0 );
    DenseVector<_Type> R_node( order, 0.0 );
    // make sure (un)refine flags vector is sized
    std::vector< bool > refine( N, false );
    std::vector< bool > unrefine( N, false );
    // Residual vector at interior nodes
    DenseVector<double> Res2( N, 0.0 );
    // reset row counter
    row = 0;
    // inner nodes of the mesh, node = 1,2,...,N-2
    for ( std::size_t node = 1; node <= N - 2; node += 2 )
    {
      // set the current solution at this node
      for ( unsigned var = 0; var < order; ++var )
      {
        //F_node[ var ] = SOLUTION[ var ][ node ];
        F_node[ var ] = SOLUTION( node, var );
      }
      // set the y-pos in the eqn
      p_EQUATION -> y() = SOLUTION.coord( node );
      // Update the equation to the nodal position
      p_EQUATION -> update( F_node );
      //// evaluate the RHS at the node
      //p_EQUATION -> get_residual( R_node );
      // step size centred at the node
      _Xtype invh = 1. / ( SOLUTION.coord( node + 1 ) - SOLUTION.coord( node - 1 ) );
      // loop over all the variables
      DenseVector<_Type> temp( order, 0.0 );
      for ( unsigned var = 0; var < order; ++var )
      {
        // residual
        //temp[ var ] = p_EQUATION -> residual()[ var ] - ( SOLUTION[ var ][ node + 1 ] - SOLUTION[ var ][ node - 1 ] ) * invh;
        temp[ var ] = p_EQUATION -> residual()[ var ] - ( SOLUTION( node + 1, var ) - SOLUTION( node - 1, var ) ) * invh;
      }
      Res2[ node ] = temp.inf_norm();
      if ( Res2[ node ] > adapt_tol )
      {
        refine[ node ] = true;
      }
      else
        if ( Res2[ node ] < TOL / 10. )
        {
          unrefine[ node ] = true;
        }
    }

    std::size_t no_refined( 0 ), no_unrefined( 0 );
    for ( std::size_t i = 0; i < refine.size(); ++i )
    {
      if ( refine[ i ] == true )
      {
        no_refined += 1;
      }
      if ( unrefine[ i ] == true )
      {
        no_unrefined += 1;
      }
    }
#ifdef DEBUG
    std::cout << "[ DEBUG ] Refinements = " << no_refined << "\n";
    std::cout << "[ DEBUG ] Unrefinements = " << no_unrefined << "\n";
#endif

    // make a new refined/unrefined mesh
    DenseVector<_Xtype> X( SOLUTION.nodes() );
    DenseVector<_Xtype> newX;
    newX.push_back( X[ 0 ] );
    for ( std::size_t i = 1; i < N - 1; ++i )
    {
      if ( refine[ i ] )
      {
        if ( !refine[ i - 1 ] )
        {
          // have not already refined to the left
          // so refine left AND right with new nodes
          _Xtype left( X[ i - 1 ] );
          _Xtype right( X[ i + 1 ] );
          _Xtype here( X[ i ] );
          newX.push_back( ( left + here ) / 2. );
          newX.push_back( here );
          newX.push_back( ( right + here ) / 2. );
        }
        else
        {
          // already left refined, so just refine right
          _Xtype right( X[ i + 1 ] );
          _Xtype here( X[ i ] );
          newX.push_back( here );
          newX.push_back( ( right + here ) / 2. );
        }
      }
      else
        if ( !unrefine[ i ] )
        {
          newX.push_back( X[ i ] );
        }
    }
    newX.push_back( X[ N - 1 ] );

    // remesh the current solution to this (un)refined mesh
    SOLUTION.remesh1( newX );
#ifdef TIME
    T_REFINE.stop();
#endif
    // adding nodes will screw up the sign of the determinant of the Jacobian
    // so here we'll just make it zero
    LAST_DET_SIGN = 0;

    std::pair< unsigned, unsigned > feedback;
    feedback.first = no_refined;
    feedback.second = no_unrefined;
    return feedback;
  }


  template <typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::assemble_linear_evp( DenseMatrix<_Type>& dense_a, DenseMatrix<_Type>& dense_b )
  {
    // cast the equation to an equation with mass
    Equation_with_mass<_Type, _Xtype>* p_masseqn = dynamic_cast< Equation_with_mass<_Type, _Xtype>* >( p_EQUATION );
    // the order of the problem
    unsigned order( p_masseqn -> get_order() );
    // get the number of nodes in the mesh
    unsigned n( SOLUTION.get_nnodes() );
    // assemble the banded/vector system for the BVP
    BandedMatrix<_Type> a( n * order, 2 * order - 1, 0.0 );
    DenseVector<_Type> b( n * order, 0.0 );
    assemble_matrix_problem( a, b );
    // move the banded matrix to a dense one
    dense_a = a;
    // make a new dense RHS - which will just be contributed to
    // by the mass matrix of the equation.
    dense_b = DenseMatrix<_Type>( n * order, n * order, 0.0 );
    // local state variable and functions
    DenseVector<_Type> F_midpt( order, 0.0 );
    // skip the boundary condition rows at the left of the domain
    std::size_t row( p_LEFT_RESIDUAL -> get_order() );
    // inner nodes of the mesh, node = 0,1,2,...,N-2
    for ( std::size_t node = 0; node <= n - 2; ++node )
    {
      const std::size_t lnode = node;
      const std::size_t rnode = node + 1;
      // set the current solution at this node by 2nd order evaluation at mid point
      for ( unsigned var = 0; var < order; ++var )
      {
        F_midpt[ var ] = ( SOLUTION( lnode, var ) + SOLUTION( rnode, var ) ) / 2.;
      }
      // mid point of the independent variable
      _Xtype y_midpt = 0.5 * ( SOLUTION.coord( lnode ) + SOLUTION.coord( rnode ) );
      // set the equation coord
      p_masseqn -> y() = y_midpt;
      // Update the equation to the mid point position
      p_masseqn -> update( F_midpt );
      // loop over all the variables
      for ( unsigned var = 0; var < order; ++var )
      {
        // add the Jacobian terms to the linearised problem
        for ( unsigned i = 0; i < order; ++i )
        {
          dense_b( row, order * lnode + i ) -= p_masseqn -> mass()( var, i ) / 2.;
          dense_b( row, order * rnode + i ) -= p_masseqn -> mass()( var, i ) / 2.;
        }
        // increment the row
        row += 1;
      }
    }
  }

  // the templated versions that we require are:
  template class ODE_BVP<double>
  ;
  template class ODE_BVP<std::complex<double> >
  ;
  template class ODE_BVP<std::complex<double>, std::complex<double> >
  ;


} // end namespace

