/// \file TwoD_TVDLF_Elt.cpp
/// Implementation of an object that represents a two dimensional
/// element for TVD LF methods.

#include <vector>

#include <list>
#include <set>
#include <algorithm>

#include <TwoD_Hyperbolic_System.h>
#include <TwoD_TVDLF_Elt.h>


namespace CppNoddy
{

  TwoD_TVDLF_Elt::TwoD_TVDLF_Elt( double west, double east,
                                  double south, double north,
                                  TwoD_Hyperbolic_System* ptr,
                                  bool flag, std::set<int> faces )
  {
#ifdef DEBUG
    std::cout << "DEBUG: Constructing a TwoD_TVDLF_Elt object over ["
              << west << ", " << east << "] X [" << south << ", " << north << "]\n";
#endif
    p_system = ptr;
    // south west corner
    this -> WEST = west;
    this -> SOUTH = south;
    // lengths of the element sides
    this -> DX = east - west;
    this -> DY = north - south;
    // is this an elt with an external face?
    EXTERNAL = flag;
    // which faces are external?
    EXTERNAL_FACES = faces;
    // pointers to the elements on the compass points
    p_ELTS = std::vector<TwoD_TVDLF_Elt*>( 4, this );
    // linear reconstruction within the elt ... in 2 directions
    SLOPE_X = DenseVector<double>( p_system -> get_order(), 0.0 );
    SLOPE_Y = DenseVector<double>( p_system -> get_order(), 0.0 );
    // the concentrations in this elt
    Q = DenseVector<double>( p_system -> get_order(), 0.0 );
    // the corrections to be added ... built up in parts
    DELTA_Q = DenseVector<double>( p_system -> get_order(), 0.0 );
    // no fluxes have been added to this elt
    FLUX_FACE_DONE = std::vector<bool>( 4, false );
  }

  TwoD_TVDLF_Elt::~TwoD_TVDLF_Elt()
  {}

  bool TwoD_TVDLF_Elt::face_is_external( const int& index ) const
  {
    bool found( true );
    if ( EXTERNAL_FACES.find( index ) == EXTERNAL_FACES.end() )
    {
      found = false;
    }
    return found;
  }

  bool TwoD_TVDLF_Elt::face_is_internal( const int& index ) const
  {
    return !face_is_external( index );
  }

  void TwoD_TVDLF_Elt::set_ptrs( const int& index, TwoD_TVDLF_Elt* ptr )
  {
    p_ELTS[ index ] = ptr;
  }

  TwoD_TVDLF_Elt* TwoD_TVDLF_Elt::get_ptrs( const int& index ) const
  {
    return p_ELTS[ index ];
  }

  void TwoD_TVDLF_Elt::add_flux_contributions( const double& dt )
  {
    DenseVector<double> flux( p_system -> get_order(), 0.0 );
    // work out the flux into the elt
    if ( !FLUX_FACE_DONE[ 0 ] )
    {
      contributed_flux_in_south( dt, flux );
    }
    if ( !FLUX_FACE_DONE[ 1 ] )
    {
      contributed_flux_out_east( dt, flux );
    }
    if ( !FLUX_FACE_DONE[ 2 ] )
    {
      contributed_flux_out_north( dt, flux );
    }
    if ( !FLUX_FACE_DONE[ 3 ] )
    {
      contributed_flux_in_west( dt, flux );
    }
    // evaluate the total integral over dt and divide by the area of the elt
    DELTA_Q *= dt / ( DX * DY );
    // include the source contribution
    DELTA_Q += get_source_contribution( dt );
    // update the concentrations
    Q += DELTA_Q;
    // reset the correction to zero
    DELTA_Q = DenseVector<double>( p_system -> get_order(), 0.0 );
    // reset the elt to having no fluxes added
    FLUX_FACE_DONE = std::vector<bool>( 4, false );
  }

  void TwoD_TVDLF_Elt::add_contribution( TwoD_TVDLF_Elt* ptr,
                                         const DenseVector<double>& sw,
                                         const DenseVector<double>& ne,
                                         std::set< int > indices )
  {
    // if any face indices are external then we ignore them
    std::set< int > internal_indices;
    set_difference( indices.begin(), indices.end(),
                    EXTERNAL_FACES.begin(), EXTERNAL_FACES.end(),
                    inserter( internal_indices, internal_indices.begin() ) );
    // loop thru the contributions to see if it is already in there
    std::list< contribution >::iterator c = CONT_LIST.begin();
    bool found = false;
    while ( c != CONT_LIST.end() )
    {
      if ( c -> elt_ptr == ptr )
      {
        found = true;
        break;
      }
      else
      {
        ++c;
      }
    }
    if ( found )
    {
      // just add face indices to the existing entry's set of faces
      c -> face_indices.insert( internal_indices.begin(), internal_indices.end() );
    }
    else
    {
      // we need to make a new contribution entry
      contribution cont;
      cont.elt_ptr = ptr;
      cont.sw = sw;
      cont.ne = ne;
      cont.face_indices.insert( internal_indices.begin(), internal_indices.end() );
      // push the contribution data into a list
      CONT_LIST.push_back( cont );
    }
  }

  void TwoD_TVDLF_Elt::dump() const
  {
    DenseVector<double> s( 2, 0.0 );
    std::cout << WEST + DX / 2 << " " << SOUTH + DY / 2 << " " << get_Q( s )[0] << "\n";
    if ( EXTERNAL_FACES.find( 1 ) != EXTERNAL_FACES.end() )
    {
      std::cout << "\n";
    }
  }

  void TwoD_TVDLF_Elt::clear_contributions()
  {
    CONT_LIST.clear();
  }

  DenseVector<double> TwoD_TVDLF_Elt::contributed_Q() const
  {
    std::list< contribution >::const_iterator c = CONT_LIST.begin();
    // start with zero
    DenseVector<double> sum( p_system -> get_order(), 0.0 );
    // loop over contributions
    while ( c != CONT_LIST.end() )
    {
      // get integral of each
      sum += c -> elt_ptr -> get_int_Q( c -> sw , c -> ne );
      ++c;
    }
    return sum;
  }

  void TwoD_TVDLF_Elt::contributed_flux_in_west( const double &dt, DenseVector<double> &flux_in_west )
  {
    DenseVector<double> flux( p_system -> get_order(), 0.0 );
    flux_in_west = DenseVector<double>( p_system -> get_order(), 0.0 );
    std::list< contribution >::iterator c = CONT_LIST.begin();
    if ( face_is_internal( 3 ) )
    {
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // true if this makes a contribution to the West face
        if ( c -> face_indices.find( 3 ) != c -> face_indices.end() )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the y-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = c -> sw[ 0 ];
          s_half[ 1 ] = ( c -> sw[1] + c -> ne[1] ) * 0.5;
          // current concentration value at the mid-time step
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_x( c -> elt_ptr -> get_x( s_half ), q_half_step, flux );
          // the length of the elements face that this computation is over (for the y-integral)
          const double sub_dy( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[1] );
          flux_in_west += flux * sub_dy;
        }
        ++c;
      }
      // flux_in_west of this elt = flux_out_east of the elt to the west
      p_ELTS[ 3 ] -> add_to_delta_Q( 1, -flux_in_west );
    }
    // if we are computing the flux thru West face and it is external
    // then we should use the user specified edge values
    else
    {
      std::list< contribution >::iterator c = CONT_LIST.begin();
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // find any contributing elts to this external face
        if ( c -> elt_ptr -> face_is_external( 3 ) )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the y-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = c -> sw[ 0 ];
          s_half[ 1 ] = ( c -> sw[1] + c -> ne[1] ) * 0.5;
          // global coordinate
          const DenseVector<double> x( c -> elt_ptr -> get_x( s_half ) );
          // get the concentration at the mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // at this point we have obtained the 'natural' boundary condition
          // but we now allow the user to overwrite this using the edge_values
          // not overwriting q_half_step means we keep the natural condition
          p_system -> edge_values( 3, x, q_half_step );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_x( x, q_half_step, flux );
          // the length of the elements face that this computation is over (for the y-integral)
          const double sub_dy( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[1] );
          flux_in_west += flux * sub_dy;
        }
        ++c;
      }
    }
    FLUX_FACE_DONE[ 3 ] = true;
    DELTA_Q += flux_in_west;
  }

  void TwoD_TVDLF_Elt::contributed_flux_in_south( const double &dt, DenseVector<double> &flux_in_south )
  {
    DenseVector<double> flux( p_system -> get_order(), 0.0 );
    flux_in_south = DenseVector<double>( p_system -> get_order(), 0.0 );
    std::list< contribution >::iterator c = CONT_LIST.begin();
    if ( face_is_internal( 0 ) )
    {
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // true if this makes a contribution to the South face
        if ( c -> face_indices.find( 0 ) != c -> face_indices.end() )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the x-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = ( c -> sw[ 0 ] + c -> ne[ 0 ] ) * 0.5;
          s_half[ 1 ] =  c -> sw[ 1 ];
          // current concentration value at the mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_y( c -> elt_ptr -> get_x( s_half ), q_half_step, flux );
          // the length of the elements face that this computation is over (for the x-integral)
          const double sub_dx( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[0] );
          flux_in_south += flux * sub_dx;
        }
        ++c;
      }
      // flux_in_south of this elt = flux_out_north of the elt to the south
      p_ELTS[ 0 ] -> add_to_delta_Q( 2, -flux_in_south );
    }
    // if we are computing the flux thru South face and it is external
    // then we should use the user specified edge values
    else
    {
      std::list< contribution >::iterator c = CONT_LIST.begin();
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // find any contributing elts to this external face
        if ( c -> elt_ptr -> face_is_external( 0 ) )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the x-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = ( c -> sw[ 0 ] + c -> ne[ 0 ] ) * 0.5;
          s_half[ 1 ] =  c -> sw[ 1 ];
          // global coordinate
          const DenseVector<double> x( c -> elt_ptr -> get_x( s_half ) );
          // get the concentration at the mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // at this point we have obtained the 'natural' boundary condition
          // but we now allow the user to overwrite this using the edge_values
          // not overwriting q_half_step means we keep the natural condition
          p_system -> edge_values( 0, x, q_half_step );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_y( x, q_half_step, flux );
          // the length of the elements face that this computation is over (for the x-integral)
          const double sub_dx( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[0] );
          flux_in_south += flux * sub_dx;
        }
        ++c;
      }
    }
    FLUX_FACE_DONE[ 0 ] = true;
    DELTA_Q += flux_in_south;
  }

  void TwoD_TVDLF_Elt::contributed_flux_out_east( const double &dt, DenseVector<double> &flux_out_east )
  {
    DenseVector<double> flux( p_system -> get_order(), 0.0 );
    flux_out_east = DenseVector<double>( p_system -> get_order(), 0.0 );
    std::list< contribution >::iterator c = CONT_LIST.begin();
    if ( face_is_internal( 1 ) )
    {
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // true if this makes a contribution to the East face
        if ( c -> face_indices.find( 1 ) != c -> face_indices.end() )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the y-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = c -> ne[ 0 ];
          s_half[ 1 ] = ( c -> sw[1] + c -> ne[1] ) * 0.5;
          // current concentration value at the mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_x( c -> elt_ptr -> get_x( s_half ), q_half_step, flux );
          // the length of the elements face that this computation is over (for the y-integral)
          const double sub_dy( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[1] );
          flux_out_east += flux * sub_dy;
        }
        ++c;
      }
      // flux_out_east of this elt = flux_in_west of the elt to the west
      p_ELTS[ 1 ] -> add_to_delta_Q( 3, flux_out_east );
    }
    // if we are computing the flux thru East face and it is external
    // then we should use the user specified edge values
    else
    {
      std::list< contribution >::iterator c = CONT_LIST.begin();
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // find any contributing elts to this external face
        if ( c -> elt_ptr -> face_is_external( 1 ) )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the y-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = c -> ne[ 0 ];
          s_half[ 1 ] = ( c -> sw[1] + c -> ne[1] ) * 0.5;
          // global coordinate
          const DenseVector<double> x( c -> elt_ptr -> get_x( s_half ) );
          // get the concentration at this point mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // at this point we have obtained the 'natural' boundary condition
          // but we now allow the user to overwrite this using the edge_values
          // not overwriting q_half_step means we keep the natural condition
          p_system -> edge_values( 1, x, q_half_step );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_x( x, q_half_step, flux );
          // the length of the elements face that this computation is over (for the y-integral)
          const double sub_dy( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[1] );
          flux_out_east += flux * sub_dy;
        }
        ++c;
      }
    }
    FLUX_FACE_DONE[ 1 ] = true;
    DELTA_Q -= flux_out_east;
  }

  void TwoD_TVDLF_Elt::contributed_flux_out_north( const double &dt, DenseVector<double> &flux_out_north )
  {
    DenseVector<double> flux( p_system -> get_order(), 0.0 );
    flux_out_north = DenseVector<double>( p_system -> get_order(), 0.0 );
    std::list< contribution >::iterator c = CONT_LIST.begin();
    if ( face_is_internal( 2 ) )
    {
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // true if this makes a contribution to the North face
        if ( c -> face_indices.find( 2 ) != c -> face_indices.end() )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the y-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = ( c -> sw[ 0 ] + c -> ne[ 0 ] ) * 0.5;
          s_half[ 1 ] =  c -> ne[ 1 ];
          // current concentration value at the mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_y( c -> elt_ptr -> get_x( s_half ), q_half_step, flux );
          // the length of the elements face that this computation is over (for the y-integral)
          const double sub_dx( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[0] );
          flux_out_north += flux * sub_dx;
        }
        ++c;
      }
      // flux_out_north of this elt = flux_in_south of the elt to the north
      p_ELTS[ 2 ] -> add_to_delta_Q( 0, flux_out_north );
    }
    // if we are computing the flux thru North face and it is external
    // then we should use the user specified edge values
    else
    {
      std::list< contribution >::iterator c = CONT_LIST.begin();
      // loop over all contributing elements
      while ( c != CONT_LIST.end() )
      {
        // find any contributing elts to this external face
        if ( c -> elt_ptr -> face_is_external( 2 ) )
        {
          // get the local coords in the contributing elt, but has to
          // be at the mid-point for the x-integration
          DenseVector<double> s_half( 2, 0.0 );
          s_half[ 0 ] = ( c -> sw[ 0 ] + c -> ne[ 0 ] ) * 0.5;
          s_half[ 1 ] =  c -> ne[ 1 ];
          // global coordinate
          const DenseVector<double> x( c -> elt_ptr -> get_x( s_half ) );
          // get the concentration at this mid-time point
          DenseVector<double> q_half_step( c -> elt_ptr -> get_Q_Taylor_series( s_half, dt / 2 ) );
          // at this point we have obtained the 'natural' boundary condition
          // but we now allow the user to overwrite this using the edge_values
          // not overwriting q_half_step means we keep the natural condition
          p_system -> edge_values( 2, x, q_half_step );
          // evaluate the flux at this Q value
          c -> elt_ptr -> p_system -> flux_fn_y( x, q_half_step, flux );
          // the length of the elements face that this computation is over (for the x-integral)
          const double sub_dx( ( c -> elt_ptr -> get_x( c -> ne ) - c -> elt_ptr -> get_x( c -> sw ) )[0] );
          flux_out_north += flux * sub_dx;
        }
        ++c;
      }
    }
    FLUX_FACE_DONE[ 2 ] = true;
    DELTA_Q -= flux_out_north;
  }

  void TwoD_TVDLF_Elt::add_to_delta_Q( const int& face_index, const DenseVector<double>& dq )
  {
    FLUX_FACE_DONE[ face_index ] = true;
    DELTA_Q += dq;
  }

  DenseVector<double> TwoD_TVDLF_Elt::get_s( const DenseVector<double> &x ) const
  {
#ifdef PARANOID
    if ( ( x[0] < WEST ) || ( x[0] > WEST + DX ) || ( x[1] > SOUTH + DY ) || ( x[1] < SOUTH ) )
    {
      std::string problem;
      problem = " The TwoD_TVDLF_Elt::get_s method has been called for an element \n";
      problem += " whose range does not bracket the given global coordinate.\n";
      throw ExceptionRuntime( problem );
    }
#endif
    DenseVector<double> s( 2, 0.0 );
    s[0] = -1 + 2 * ( x[0] - WEST ) / DX;
    s[1] = -1 + 2 * ( x[1] - SOUTH ) / DY;
    return s;
  }

  DenseVector<double> TwoD_TVDLF_Elt::get_x( DenseVector<double> s ) const
  {
    DenseVector<double> x( 2, 0.0 );
    x[ 0 ] = WEST + ( s[ 0 ] + 1 ) * DX / 2;
    x[ 1 ] = SOUTH + ( s[ 1 ] + 1 ) * DY / 2;
    return x;
  }

  DenseVector<double> TwoD_TVDLF_Elt::get_x_mid() const
  {
    DenseVector<double> x( 2, 0.0 );
    x[ 0 ] = WEST + DX / 2;
    x[ 1 ] = SOUTH + DY / 2;
    return x;
  }

  DenseVector<double> TwoD_TVDLF_Elt::get_dx() const
  {
    DenseVector<double> delta( 2, 0.0 );
    delta[ 0 ] = DX;
    delta[ 1 ] = DY;
    return delta;
  }

  double TwoD_TVDLF_Elt::get_dA() const
  {
    return DX * DY;
  }

  void TwoD_TVDLF_Elt::set_external_flag( bool flag )
  {
    EXTERNAL = flag;
  }

  bool TwoD_TVDLF_Elt::get_external_flag() const
  {
    return EXTERNAL;
  }

  std::set<int> TwoD_TVDLF_Elt::get_external_faces() const
  {
    return EXTERNAL_FACES;
  }

  void TwoD_TVDLF_Elt::add_external_face( const int& i )
  {
    EXTERNAL_FACES.insert( i );
  }

  void TwoD_TVDLF_Elt::remove_external_face( const int& i )
  {
    EXTERNAL_FACES.erase( i );
  }

  DenseVector<double> TwoD_TVDLF_Elt::get_Q_Taylor_series( const DenseVector<double> &s, const double &dt ) const
  {
    const DenseVector<double> x( get_x( s ) );
    // value of Q at the current time point
    DenseVector<double> q( get_Q( s ) );
    // evaluate the 2nd order Taylor series expansion of Q at the
    // time step of size dt
    //
    // get the RHS from the eqn
    DenseVector<double> r( p_system -> get_order(), 0.0 );
    p_system -> source_fn( x, q, r );
    // get the Jacobian of the x-flux fn
    DenseMatrix<double> jac( p_system -> get_order(), p_system -> get_order(), 0.0 );
    p_system -> Jac_flux_fn_x( x, q, jac );
    r -= jac.multiply( SLOPE_X );
    // get the Jacobian of the y-flux fn
    p_system -> Jac_flux_fn_y( x, q, jac );
    r -= jac.multiply( SLOPE_Y );
    // add the correction
    q += r * dt;
    // the above should be equivalent to the expression below, but
    // minimises instantiation of Jacobian matrix
    // q += ( get_source_fn( s )
    //   - get_Jac_flux_fn_x( s ).multiply( slope_x )
    //   - get_Jac_flux_fn_y( s ).multiply( slope_y ) ) * dt;
    return q;
  }

  DenseVector<double> TwoD_TVDLF_Elt::get_source_contribution( const double &dt ) const
  {
    const DenseVector<double> s_mid( 2, 0.0 );
    // evaluate the Taylor series expansion of the unknowns at the
    // mid time displacement point dt/2
    DenseVector<double> q( get_Q_Taylor_series( s_mid, dt / 2 ) );
    DenseVector<double> R( p_system -> get_order(), 0.0 );
    // work out the source function for this mid-point in time q value
    p_system -> source_fn( get_x( s_mid ), q, R );
    // return the contribution
    return R * dt;
  }



}
