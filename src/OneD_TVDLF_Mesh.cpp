/// \file OneD_TVDLF_Mesh.cpp
/// Implementation of an object that represents a one dimensional
/// mesh for TVD LF methods.

#include <vector>

#include <Types.h>
#include <OneD_TVDLF_Elt.h>
#include <OneD_Node_Mesh.h>
#include <OneD_Hyperbolic_System.h>
#include <OneD_TVDLF_Mesh.h>

namespace CppNoddy
{
  /// Constructor for the Finite Volume Mesh using linear elements
  /// \param X A vector of nodal locations at which the element
  ///   FACES will positioned
  OneD_TVDLF_Mesh::OneD_TVDLF_Mesh( const DenseVector<double>& X,
                                    OneD_Hyperbolic_System* ptr,
                                    fn_ptr init_ptr )
  {
    MESH_TIME = 0.0;
#ifdef DEBUG
    std::cout << "DEBUG: Starting construction of a OneD_TVDLF_Mesh object. \n";
#endif
    unsigned N = X.size();
    if ( N <= 2 )
    {
      std::string problem;
      problem = " The OneD_TVDLF_Mesh object is trying to construct itself \n";
      problem += " with just one element! \n";
      throw ExceptionRuntime( problem );
    }
#ifdef DEBUG
    std::cout << "\nDEBUG: configuration of the black mesh \n";
#endif
    // set up the fn ptr to the initial conditions fn
    p_Q_INIT = init_ptr;
    // store the order of the conservative system here for simplicity
    ORDER_OF_SYSTEM = ptr -> get_order();

    // set up the black elements
    {
      // first elt
      BLACK_ELTS.push_back( OneD_TVDLF_Elt( X[0], X[1], ptr, true, -1 ) );
      for ( std::size_t i = 1; i <= N - 3; ++i )
      {
        // interior elts
        BLACK_ELTS.push_back( OneD_TVDLF_Elt( X[i], X[i+1], ptr ) );
      }
      // last elt
      BLACK_ELTS.push_back( OneD_TVDLF_Elt( X[N-2], X[N-1], ptr, true, 1 ) );
    }
#ifdef DEBUG
    std::cout << "\nDEBUG: configuration of the red mesh \n";
#endif
    // set up the red elements
    {
      // first elt
      RED_ELTS.push_back( OneD_TVDLF_Elt( X[0], ( X[1] + X[0] ) / 2, ptr, true, -1 ) );
      // interior elts
      for ( std::size_t i = 1; i <= N - 2; ++i )
      {
        // interior elts
        RED_ELTS.push_back( OneD_TVDLF_Elt( ( X[i-1] + X[i] ) / 2, ( X[i] + X[i+1] ) / 2, ptr ) );
      }
      // last elt
      RED_ELTS.push_back( OneD_TVDLF_Elt( ( X[N-2] + X[N-1] ) / 2, X[N-1], ptr, true, 1 ) );
    }
#ifdef DEBUG
    std::cout << "\nDEBUG: linking the two meshes \n";
#endif
    // the only tricky part for this scheme is the mesh interconnections
    // black to red is easy enough
    elt_iter er = RED_ELTS.begin();
    elt_iter eb = BLACK_ELTS.begin();
    DenseVector<double> s_left( 2, 0. );
    s_left[ 0 ] = -1.;
    s_left[ 1 ] = 0.;
    DenseVector<double> s_right( 2, 0. );
    s_right[ 0 ] = 0.;
    s_right[ 1 ] = 1.;
    DenseVector<double> s_whole( 2, 0. );
    s_whole[ 0 ] = -1.;
    s_whole[ 1 ] = 1.;
    DenseVector<double> s_gen( 2, 0. );
    // loop over red elts -- and define which black elts contribute
    // this is straightforward, even for nonuniform meshes.
    while ( er != RED_ELTS.end() )
    {
      if ( er -> get_external_flag() )
      {
        // if its an external elt
        if ( er -> get_external_face_i() < 0 )
        {
          // and external face is left
          er -> add_contribution( &( *eb ), s_left, 1 );
          ++er;
        }
        else
        {
          // and external face is right
          er -> add_contribution( &( *eb ), s_right, -1 );
          ++er;
        }
      }
      else
      {
        //internal elt
        er -> add_contribution( &( *eb ), s_right, -1 );
        ++eb;
        er -> add_contribution( &( *eb ), s_left, 1 );
        ++er;
      }
    }
    // loop over all black elts and define which red elts contribute
    // this is more involved if we allow for non-uniform meshes.
    eb = BLACK_ELTS.begin();
    er = RED_ELTS.begin();
    while ( eb != BLACK_ELTS.end() )
    {
      if ( eb -> get_external_flag() )
      {
        // if its an external elt
        if ( eb -> get_external_face_i() < 0 )
        {
          // and external face is left
          eb -> add_contribution( &( *er ), s_whole, 0 );
          ++er;
          double s = er -> get_s( eb -> get_x( 1.0 ) );
          s_gen[ 0 ] = -1.0;
          s_gen[ 1 ] = s;
          eb -> add_contribution( &( *er ), s_gen, 1 );
          ++eb;
        }
        else
        {
          // and external face is right
          double s = er -> get_s( eb -> get_x( -1.0 ) );
          s_gen[ 0 ] = s;
          s_gen[ 1 ] = 1.0;
          eb -> add_contribution( &( *er ), s_gen, -1 );
          ++er;
          eb -> add_contribution( &( *er ), s_whole, 0 );
          ++eb;
        }
      }
      else
      {
        // internal elt
        double s = er -> get_s( eb -> get_x( -1.0 ) );
        s_gen[ 0 ] = s;
        s_gen[ 1 ] = 1.0;
        eb -> add_contribution( &( *er ), s_gen, -1 );
        ++er;
        s = er -> get_s( eb -> get_x( 1.0 ) );
        s_gen[ 0 ] = -1.0;
        s_gen[ 1 ] = s;
        eb -> add_contribution( &( *er ), s_gen, 1 );
        ++eb;
      }
    }

#ifdef DEBUG
    std::cout << "DEBUG: Setting the initial state of the meesh. \n";
#endif
    // set the initial condition for each elt
    eb = BLACK_ELTS.begin();
    while ( eb != BLACK_ELTS.end() )
    {
      DenseVector<double> Q( ORDER_OF_SYSTEM, 0.0 );
      p_Q_INIT( eb -> get_x( 0.0 ), Q );
      eb -> set_Q_mid( Q );
      ++eb;
    }

    //eb = BLACK_ELTS.end();
    //--eb;
    //std::cout << "Last elt = " << eb -> get_Q( 0.0 )[ 1 ] << "\n";
    //std::cout << "size = " << eb -> get_dx() << "\n";
    //--eb;
    //std::cout << "Last elt = " << eb -> get_Q( 0.0 )[ 1 ] << "\n";
    //std::cout << "size = " << eb -> get_dx() << "\n";

    //DenseVector<double> x( get_mid_node_vector() );
    //OneD_Mesh<double> soln( x );
    //soln.set_nvars( ORDER_OF_SYSTEM );
    //for ( std::size_t i = 0; i < x.size(); ++i )
    //{
    //soln.set_nodes_vars( i, BLACK_ELTS[i].get_Q( 0.0 ) );
    //std::cout << "Get_soln Q = " << soln.get_nodes_vars( i )[ 1 ] << " at x = " << x[i] << "\n";
    //}

    // default limiter = 0 (minmod)
    LIMITER = 0;
    calc_slopes( &BLACK_ELTS );

#ifdef DEBUG
    std::cout << "DEBUG: Finished construction of a OneD_TVDLF_Mesh object. \n";
#endif
  }

  OneD_TVDLF_Mesh::~OneD_TVDLF_Mesh()
  {}

  DenseVector<double> OneD_TVDLF_Mesh::get_mid_node_vector() const
  {
    DenseVector<double> X;
    celt_iter e = BLACK_ELTS.begin();
    while ( e != BLACK_ELTS.end() )
    {
      X.push_back( e -> get_x( 0.0 ) );
      ++e;
    }
    return X;
  }

  DenseVector<double> OneD_TVDLF_Mesh::get_face_pos_vector() const
  {
    DenseVector<double> X;
    celt_iter e = BLACK_ELTS.begin();
    while ( e != BLACK_ELTS.end() )
    {
      X.push_back( e -> get_x( -1.0 ) );
      ++e;
    }
    --e;
    X.push_back( e -> get_x( 1.0 ) );
    return X;
  }

  void OneD_TVDLF_Mesh::set_limiter( unsigned id )
  {
    LIMITER = id;
  }

  double OneD_TVDLF_Mesh::update( const double& CFL, const double& max_dt )
  {
    double first_dt = update_to_red( CFL, max_dt / 2.0 );
    double second_dt = update_to_black( CFL, max_dt - first_dt );
    return first_dt + second_dt;
  }

  void OneD_TVDLF_Mesh::update_to( const double& CFL, const double& t_end )
  {
    do
    {
      update( CFL, std::abs( t_end - MESH_TIME ) );
    }
    while ( MESH_TIME < t_end );
  }

  const double& OneD_TVDLF_Mesh::get_time() const
  {
    return MESH_TIME;
  }

  DenseVector<double> OneD_TVDLF_Mesh::integrate() const
  {
    celt_iter e = BLACK_ELTS.begin();
    DenseVector<double> sum( ORDER_OF_SYSTEM, 0.0 );
    while ( e != BLACK_ELTS.end() )
    {
      sum += e -> get_Q( 0.0 ) * e -> get_dx();
      ++e;
    }
    return sum;
  }

  OneD_Node_Mesh<double> OneD_TVDLF_Mesh::get_soln( std::string location, std::string mesh_colour )
  {
    vector_of_elts* elts( get_elts_from_colour( mesh_colour ) );
    OneD_Node_Mesh<double> soln;
    if ( location == "centre" )
    {
      DenseVector<double> X( get_mid_node_vector() );
      soln = OneD_Node_Mesh<double>( X, ORDER_OF_SYSTEM );
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        soln.set_nodes_vars( i, ( *elts )[i].get_Q( 0.0 ) );
      }
    }
    else
      if ( location == "face_average" )
      {
        DenseVector<double> X( get_face_pos_vector() );
        soln = OneD_Node_Mesh<double>( X, ORDER_OF_SYSTEM );
        std::size_t N = X.size() - 1;
        soln.set_nodes_vars( 0, ( *elts )[0].get_Q( -1.0 ) );
        for ( std::size_t i = 1; i < N; ++i )
        {
          soln.set_nodes_vars( i, ( ( *elts )[i-1].get_Q( 1.0 ) + ( *elts )[i].get_Q( -1.0 ) ) / 2 );
        }
        soln.set_nodes_vars( N, ( *elts )[N-1].get_Q( 1.0 ) );
      }
      else
      {
        std::string problem;
        problem = " In OneD_TVDLF_Mesh::get_soln you have passed an unrecognised ";
        problem += " location for data output. Use 'centre' or 'face_average'. \n";
        throw ExceptionRuntime( problem );
      }

    return soln;
  }

  void OneD_TVDLF_Mesh::calc_slopes( vector_of_elts* elt_vector )
  {
    set_boundary_Q( elt_vector );
    elt_iter e = elt_vector -> begin();
    DenseVector<double> slope( ORDER_OF_SYSTEM, 0.0 );
    // having computed the right-slope in elt-i we know the left-slope in elt-i+1
    DenseVector<double> left = left_diff( e );
    DenseVector<double> right( ORDER_OF_SYSTEM, 0.0 );
    switch ( LIMITER )
    {
    case 0:
      while ( e != elt_vector -> end() )
      {
        // minmod
        right = right_diff( e );
        slope = minmod( left, right );
        e -> set_slope( slope );
        left = right;
        ++e;
      }
      break;
    case 1:
      while ( e != elt_vector -> end() )
      {
        right = right_diff( e );
        slope = maxmod(
                  minmod( right, left * 2. ),
                  minmod( right * 2., left )
                );
        left = right;
        e -> set_slope( slope );
        ++e;
      }
      break;
    default:
      std::string problem;
      problem = " The OneD_TVDLF_Mesh object has an unrecognised 'limiter' identifier. \n";
      throw ExceptionRuntime( problem );
    }
  }

  std::vector<OneD_TVDLF_Elt>*  OneD_TVDLF_Mesh::get_elts_from_colour( std::string mesh_colour )
  {
    vector_of_elts* elts;
    if ( mesh_colour == "black" )
    {
      elts = &BLACK_ELTS;
    }
    else
      if ( mesh_colour == "red" )
      {
        elts = &RED_ELTS;
      }
      else
      {
        std::string problem;
        problem = " The TwoD_TVDLF_Mesh object has been passed an unrecognised mesh identifier. \n";
        problem += " valid options are 'black' or 'red'.\n";
        throw ExceptionRuntime( problem );
      }
    return elts;
  }

  void OneD_TVDLF_Mesh::set_boundary_Q( vector_of_elts* elts )
  {
    // loop over all elts
    elt_iter e( elts -> begin() );
    while ( e != elts -> end() )
    {
      // only consider the external elts
      if ( e -> get_external_flag() )
      {
        // get local coord and index of the external 'face'
        double s;
        int face;
        if ( e -> get_external_face_i() < 0 )
        {
          // left is external
          face = -1;
          s = -1.0;
        }
        else
        {
          // right is external
          face = 1;
          s = 1.0;
        }
        // get the edge value
        DenseVector<double> Qe( e -> get_Q( s ) );
        // allow the user to overwrite this value
        std::vector<bool> inflow = e -> system_ptr -> edge_values( face, e -> get_x( s ) , Qe );
        DenseVector<double> Qm( e -> get_Q( 0.0 ) );
        // if the value has been specified as inflow
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          // make the nodal value the edge value
          if ( inflow[ i ] == true )
          {
            Qm[ i ] = Qe[ i ];
          }
        }
        e -> set_Q_mid( Qm );
      }
      ++e;
    }
  }

  DenseVector<double> OneD_TVDLF_Mesh::left_diff( elt_iter e )
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    if ( e -> get_external_face_i() < 0 )
    {
      // interior difference
      diff = right_diff( e );
      // now set the deriv to zero for any inflow BCs
      //
      // get the edge value
      DenseVector<double> Qe( e -> get_Q( -1.0 ) );
      // allow the user to overwrite this value
      std::vector<bool> inflow = e -> system_ptr -> edge_values( -1, e -> get_x( -1.0 ), Qe );
      // if the value has been specified as inflow
      for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
      {
        // set a zero slope as per Levy & Tadmor (1997)
        if ( inflow[ i ] == true )
        {
          diff[ i ] = 0.0;
        }
      }
    }
    else
    {
      // We're not on an edge, so we can difference
      elt_iter j( e );
      --j;
      diff = ( e -> get_Q( 0.0 ) - j -> get_Q( 0.0 ) )
             / ( e -> get_x( 0.0 ) - j -> get_x( 0.0 ) );
    }
    return diff;
  }

  DenseVector<double> OneD_TVDLF_Mesh::right_diff( elt_iter e )
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    //if ( ( e -> get_external_flag() ) && ( e -> get_external_face_i() > 0 ) )
    if ( e -> get_external_face_i() > 0 )
    {
      // interior difference
      diff = left_diff( e );
      // now set the deriv to zero for any inflow BCs
      //
      // get the edge value
      DenseVector<double> Qe( e -> get_Q( 1.0 ) );
      // allow the user to overwrite this value
      std::vector<bool> inflow = e -> system_ptr -> edge_values( 1, e -> get_x( 1.0 ), Qe );
      DenseVector<double> Qm( e -> get_Q( 0.0 ) );
      // if the value has been specified as inflow
      for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
      {
        // set a zero slope as per Levy & Tadmor (1997)
        if ( inflow[ i ] == true )
        {
          diff[ i ] = 0.0;
        }
      }
    }
    else
    {
      // We're not on an edge, so we can difference
      elt_iter j( e );
      ++j;
      diff = ( j -> get_Q( 0.0 ) - e -> get_Q( 0.0 ) )
             / ( j -> get_x( 0.0 ) - e -> get_x( 0.0 ) );
    }
    return diff;
  }

  int OneD_TVDLF_Mesh::sgn( double a )
  {
    if ( a > 0.0 )
    {
      return 1;
    }
    else
      if ( a < 0.0 )
      {
        return -1;
      }
      else
      {
        return 0;
      }
  }

  double OneD_TVDLF_Mesh::minmod( double a, double b )
  {
    return 0.5 * ( sgn( a ) + sgn( b ) )
           * std::min( std::abs( a ), std::abs( b ) );
  }

  double OneD_TVDLF_Mesh::maxmod( double a, double b )
  {
    return 0.5 * ( sgn( a ) + sgn( b ) )
           * std::max( std::abs( a ), std::abs( b ) );
  }

  DenseVector<double> OneD_TVDLF_Mesh::minmod( DenseVector<double> A, DenseVector<double> B )
  {
    DenseVector<double> MM;
    for ( std::size_t i = 0; i < A.size(); ++i )
    {
      MM.push_back( 0.5 * ( sgn( A[i] ) + sgn( B[i] ) )
                    * std::min( std::abs( A[i] ), std::abs( B[i] ) ) );
    }
    return MM;
  }

  DenseVector<double> OneD_TVDLF_Mesh::maxmod( DenseVector<double> A, DenseVector<double> B )
  {
    DenseVector<double> MM;
    for ( std::size_t i = 0; i < A.size(); ++i )
    {
      MM.push_back( 0.5 * ( sgn( A[i] ) + sgn( B[i] ) )
                    * std::max( std::abs( A[i] ), std::abs( B[i] ) ) );
    }
    return MM;
  }

}
