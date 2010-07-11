/// \file TwoD_TVDLF_Mesh.cpp
/// Implementation of an object that represents a two dimensional
/// mesh for TVD LF methods.

#include <vector>
#include <algorithm>

#include <Types.h>
#include <TwoD_TVDLF_Elt.h>
#include <TwoD_Hyperbolic_System.h>
#include <TwoD_TVDLF_Mesh.h>

namespace CppNoddy
{

  TwoD_TVDLF_Mesh::TwoD_TVDLF_Mesh( const DenseVector<double>& X, const DenseVector<double>& Y,
                                    TwoD_Hyperbolic_System* ptr,
                                    fn_ptr init_ptr )
  {

#ifdef DEBUG
    std::cout << "DEBUG: Starting construction of a TwoD_TVDLF_Mesh object. \n";
#endif
    MESH_TIME = 0.0;
    NX = X.size();
    NY = Y.size();
    if ( std::min( NX - 1, NY - 1 ) <= 1 )
    {
      std::string problem;
      problem = " The TwoD_TVDLF_Mesh object is trying to construct itself \n";
      problem += " with just one element in one of the directions! \n";
      throw ExceptionRuntime( problem );
    }

#ifdef DEBUG
    std::cout << "DEBUG: configuration of the black mesh \n";
#endif

    // reserve appropriate space to avoid re-allocation later
    BLACK_ELTS.reserve( ( NX - 1 ) * ( NY - 1 ) );
    RED_ELTS.reserve( NX * NY );

    // set up the fn ptr to the initial conditions fn
    p_Q_INIT = init_ptr;
    // default limiter
    LIMITER = 0;
    // store the order of the conservative system here for simplicity
    ORDER_OF_SYSTEM = ptr -> get_order();

    // face index of edges
    std::set<int> faces_s;
    faces_s.insert( 0 );
    std::set<int> faces_e;
    faces_e.insert( 1 );
    std::set<int> faces_n;
    faces_n.insert( 2 );
    std::set<int> faces_w;
    faces_w.insert( 3 );
    // face indices of corners
    std::set<int> faces_se( faces_s );
    faces_se.insert( 1 );
    std::set<int> faces_ne( faces_n );
    faces_ne.insert( 1 );
    std::set<int> faces_nw( faces_n );
    faces_nw.insert( 3 );
    std::set<int> faces_sw( faces_s );
    faces_sw.insert( 3 );
    // set up the black elements
    {
      // southern row
      BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[0], X[1], Y[0], Y[1], ptr, true, faces_sw ) );
      for ( std::size_t i = 1; i <= NX - 3; ++i )
      {
        BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[i], X[i+1], Y[0], Y[1], ptr, true, faces_s ) );
      }
      BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[NX-2], X[NX-1], Y[0], Y[1], ptr, true, faces_se ) );
      // interior elts
      for ( std::size_t j = 1; j <= NY - 3; ++j )
      {
        BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[0], X[1], Y[j], Y[j+1], ptr, true, faces_w ) );
        for ( std::size_t i = 1; i <= NX - 3; ++i )
        {
          BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[i], X[i+1], Y[j], Y[j+1], ptr ) );
        }
        BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[NX-2], X[NX-1], Y[j], Y[j+1], ptr, true, faces_e ) );
      }
      // northern row
      BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[0], X[1], Y[NY-2], Y[NY-1], ptr, true, faces_nw ) );
      for ( std::size_t i = 1; i <= NX - 3; ++i )
      {
        BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[i], X[i+1], Y[NY-2], Y[NY-1], ptr, true, faces_n ) );
      }
      BLACK_ELTS.push_back( TwoD_TVDLF_Elt( X[NX-2], X[NX-1], Y[NY-2], Y[NY-1], ptr, true, faces_ne ) );
    }

#ifdef DEBUG
    std::cout << "DEBUG: configuration of the red mesh \n";
#endif

    // set up the red elements
    {
      // southern row
      RED_ELTS.push_back( TwoD_TVDLF_Elt( X[0], ( X[1] + X[0] ) / 2, Y[0], ( Y[0] + Y[1] ) / 2, ptr, true, faces_sw ) );
      for ( std::size_t i = 1; i <= NX - 2; ++i )
      {
        RED_ELTS.push_back( TwoD_TVDLF_Elt( ( X[i-1] + X[i] ) / 2, ( X[i] + X[i+1] ) / 2, Y[0], ( Y[0] + Y[1] ) / 2, ptr, true, faces_s  ) );
      }
      RED_ELTS.push_back( TwoD_TVDLF_Elt( ( X[NX-2] + X[NX-1] ) / 2, X[NX-1], Y[0], ( Y[0] + Y[1] ) / 2, ptr, true, faces_se ) );
      // interior elts
      for ( std::size_t j = 1; j <= NY - 2; ++j )
      {
        RED_ELTS.push_back( TwoD_TVDLF_Elt( X[0], ( X[1] + X[0] ) / 2, ( Y[j-1] + Y[j] ) / 2, ( Y[j] + Y[j+1] ) / 2, ptr, true, faces_w ) );
        for ( std::size_t i = 1; i <= NX - 2; ++i )
        {
          RED_ELTS.push_back( TwoD_TVDLF_Elt( ( X[i-1] + X[i] ) / 2, ( X[i] + X[i+1] ) / 2, ( Y[j-1] + Y[j] ) / 2, ( Y[j] + Y[j+1] ) / 2, ptr ) );
        }
        RED_ELTS.push_back( TwoD_TVDLF_Elt( ( X[NX-2] + X[NX-1] ) / 2, X[NX-1], ( Y[j-1] + Y[j] ) / 2, ( Y[j] + Y[j+1] ) / 2, ptr, true, faces_e ) );
      }
      // northern row
      RED_ELTS.push_back( TwoD_TVDLF_Elt( X[0], ( X[1] + X[0] ) / 2, ( Y[NY-2] + Y[NY-1] ) / 2, Y[NY-1], ptr, true, faces_nw ) );
      for ( std::size_t i = 1; i <= NX - 2; ++i )
      {
        RED_ELTS.push_back( TwoD_TVDLF_Elt( ( X[i-1] + X[i] ) / 2, ( X[i] + X[i+1] ) / 2, ( Y[NY-2] + Y[NY-1] ) / 2, Y[NY-1], ptr, true, faces_n ) );
      }
      RED_ELTS.push_back( TwoD_TVDLF_Elt( ( X[NX-2] + X[NX-1] ) / 2, X[NX-1], ( Y[NY-2] + Y[NY-1] ) / 2, Y[NY-1], ptr, true, faces_ne ) );
    }

#ifdef DEBUG
    std::cout << "DEBUG: computing black to red projection \n";
#endif

    // Now we need to pre-compute the mesh projection from black-red-black.
    // Each element needs to know which element in the other mesh contributes
    // to it in the projection scheme. It also needs to know which of its faces
    // the element contributes to in the flux computation.
    // We assume the constructor is not required to be efficient & just brute
    // force it here rather than working out all the index mappings by hand.
    //
    // loop through all the black elts
    elt_iter eb = BLACK_ELTS.begin();
    while ( eb != BLACK_ELTS.end() )
    {
      // local coordinates of the corners
      DenseVector<double> ne( 2, 1.0 );
      DenseVector<double> sw( -ne );
      DenseVector<double> se( 2, 1.0 );
      se[ 1 ] = -1.0;
      DenseVector<double> nw( -se );
      // global coordinates of the corners
      DenseVector<double> bx_sw = eb -> get_x( sw );
      DenseVector<double> bx_ne = eb -> get_x( ne );
      DenseVector<double> bx_nw = eb -> get_x( nw );
      DenseVector<double> bx_se = eb -> get_x( se );
      // get an iterator to the red elt that contains each corner
      // & add it as a contribution
      //
      // nw corner
      {
        //elt_iter er = get_elt_iter_from_x( bx_nw, "red" );
        elt_iter er = get_elt_iter_from_elt_iter( eb, "red", 3 );
        DenseVector<double> s = er -> get_s( bx_nw );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = s[ 0 ];
        s1[ 1 ] = -1.0;
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = 1.0;
        s2[ 1 ] = s[ 1 ];
        eb -> add_contribution( &( *er ), s1, s2, faces_nw );
      }
      // sw corner
      {
        //elt_iter er = get_elt_iter_from_x( bx_sw, "red" );
        elt_iter er = get_elt_iter_from_elt_iter( eb, "red", 0 );
        DenseVector<double> s = er -> get_s( bx_sw );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = s[ 0 ];
        s1[ 1 ] = s[ 1 ];
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = 1.0;
        s2[ 1 ] = 1.0;
        eb -> add_contribution( &( *er ), s1, s2, faces_sw );
      }
      // ne corner
      {
        //elt_iter er = get_elt_iter_from_x( bx_ne, "red" );
        elt_iter er = get_elt_iter_from_elt_iter( eb, "red", 2 );
        DenseVector<double> s = er -> get_s( bx_ne );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = -1.0;
        s1[ 1 ] = -1.0;
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = s[ 0 ];
        s2[ 1 ] = s[ 1 ];
        eb -> add_contribution( &( *er ), s1, s2, faces_ne );
      }
      // se corner
      {
        //elt_iter er = get_elt_iter_from_x( bx_se, "red" );
        elt_iter er = get_elt_iter_from_elt_iter( eb, "red", 1 );
        DenseVector<double> s = er -> get_s( bx_se );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = -1.0;
        s1[ 1 ] = s[ 1 ];
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = s[ 0 ];
        s2[ 1 ] = 1.0;
        eb -> add_contribution( &( *er ), s1, s2, faces_se );
      }
      ++eb;
    }

#ifdef DEBUG
    std::cout << "DEBUG: computing red to black projection \n";
#endif

    // now we need to set up the projection from red to black elts.
    // it's a little more tricky here as we have to avoid the external nodes.
    elt_iter er = RED_ELTS.begin();
    while ( er != RED_ELTS.end() )
    {
      // local coordinates of the corners
      DenseVector<double> ne( 2, 1.0 );
      DenseVector<double> sw( -ne );
      DenseVector<double> se( 2, 1.0 );
      se[ 1 ] = -1.0;
      DenseVector<double> nw( -se );
      // global coordinates of the corners
      DenseVector<double> bx_sw = er -> get_x( sw );
      DenseVector<double> bx_ne = er -> get_x( ne );
      DenseVector<double> bx_nw = er -> get_x( nw );
      DenseVector<double> bx_se = er -> get_x( se );
      //
      std::set<int> external = er -> get_external_faces();
      if ( ( external.find( 2 ) == external.end() ) &&
           ( external.find( 3 ) == external.end() ) )
      {
        // nw corner is internal
        //elt_iter eb = get_elt_iter_from_x( bx_nw, "black" );
        elt_iter eb = get_elt_iter_from_elt_iter( er, "black", 3 );
        DenseVector<double> s = eb -> get_s( bx_nw );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = s[ 0 ];
        s1[ 1 ] = -1.0;
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = 1.0;
        s2[ 1 ] = s[ 1 ];
        er -> add_contribution( &( *eb ), s1, s2, faces_nw );
      }
      if ( ( external.find( 0 ) == external.end() ) &&
           ( external.find( 3 ) == external.end() ) )
      {
        // sw corner is internal
        //elt_iter eb = get_elt_iter_from_x( bx_sw, "black" );
        elt_iter eb = get_elt_iter_from_elt_iter( er, "black", 0 );
        DenseVector<double> s = eb -> get_s( bx_sw );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = s[ 0 ];
        s1[ 1 ] = s[ 1 ];
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = 1.0;
        s2[ 1 ] = 1.0;
        er -> add_contribution( &( *eb ), s1, s2, faces_sw );
      }
      if ( ( external.find( 1 ) == external.end() ) &&
           ( external.find( 2 ) == external.end() ) )
      {
        // ne corner is internal
        //elt_iter eb = get_elt_iter_from_x( bx_ne, "black" );
        elt_iter eb = get_elt_iter_from_elt_iter( er, "black", 2 );
        DenseVector<double> s = eb -> get_s( bx_ne );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = -1.0;
        s1[ 1 ] = -1.0;
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = s[ 0 ];
        s2[ 1 ] = s[ 1 ];
        er -> add_contribution( &( *eb ), s1, s2, faces_ne );
      }
      if ( ( external.find( 0 ) == external.end() ) &&
           ( external.find( 1 ) == external.end() ) )
      {
        // se corner is internal
        //elt_iter eb = get_elt_iter_from_x( bx_se, "black" );
        elt_iter eb = get_elt_iter_from_elt_iter( er, "black", 1 );
        DenseVector<double> s = eb -> get_s( bx_se );
        DenseVector<double> s1( 2, 0.0 );
        s1[ 0 ] = -1.0;
        s1[ 1 ] = s[ 1 ];
        DenseVector<double> s2( 2, 0.0 );
        s2[ 0 ] = s[ 0 ];
        s2[ 1 ] = 1.0;
        er -> add_contribution( &( *eb ), s1, s2, faces_se );
      }
      ++er;
    }

#ifdef DEBUG
    std::cout << "DEBUG: setting up pointers to NESW elts in black mesh\n";
#endif

    // we need to set the pointers to neighbouring elts in both meshes
    eb = BLACK_ELTS.begin();
    while ( eb != BLACK_ELTS.end() )
    {
      // the offset for computing N & S elts is no. of faces minus one
      const std::size_t offset( NX - 1 );
      // get a set of external faces
      std::set< int > faces( eb -> get_external_faces() );
      if ( eb -> face_is_internal( 0 ) )
      {
        // southern face is internal
        eb -> set_ptrs( 0, &( *eb ) - offset );
      }
      if ( eb -> face_is_internal( 1 ) )
      {
        // eastern face is internal
        eb -> set_ptrs( 1, &( *eb ) + 1 );
      }
      if ( eb -> face_is_internal( 2 ) )
      {
        // northern face is internal
        eb -> set_ptrs( 2, &( *eb ) + offset );
      }
      if ( eb -> face_is_internal( 3 ) )
      {
        // western face is internal
        eb -> set_ptrs( 3, &( *eb ) - 1 );
      }
      ++eb;
    }

#ifdef DEBUG
    std::cout << "DEBUG: setting up pointers to NESW elts in red mesh\n";
#endif

    // we need to set the pointers to neighbouring elts in both meshes
    er = RED_ELTS.begin();
    while ( er != RED_ELTS.end() )
    {
      // the offset for computing N & S elts is no. of faces
      // as the red mesh has an extra elt per row
      const std::size_t offset( NX );
      // get a set of external faces
      std::set< int > faces( er -> get_external_faces() );
      if ( er -> face_is_internal( 0 ) )
      {
        // southern face is internal
        er -> set_ptrs( 0, &( *er ) - offset );
      }
      if ( er -> face_is_internal( 1 ) )
      {
        // eastern face is internal
        er -> set_ptrs( 1, &( *er ) + 1 );
      }
      if ( er -> face_is_internal( 2 ) )
      {
        // northern face is internal
        er -> set_ptrs( 2, &( *er ) + offset );
      }
      if ( er -> face_is_internal( 3 ) )
      {
        // western face is internal
        er -> set_ptrs( 3, &( *er ) - 1 );
      }
      ++er;
    }

#ifdef DEBUG
    std::cout << "DEBUG: initialising the black mesh \n";
#endif

    // now all that remains is to initialise the starting (black) mesh
    eb = BLACK_ELTS.begin();
    while ( eb != BLACK_ELTS.end() )
    {
      DenseVector<double> s( 2, 0.0 );
      // compute to the east
      s[ 0 ] = 1.0;
      s[ 1 ] = 0.0;
      DenseVector<double> xe( eb -> get_x( s ) );
      DenseVector<double> Qe( ORDER_OF_SYSTEM, 0.0 );
      p_Q_INIT( xe[0], xe[1], Qe );
      // compute to the west
      s[ 0 ] = -1.0;
      s[ 1 ] = 0.0;
      DenseVector<double> xw( eb -> get_x( s ) );
      DenseVector<double> Qw( ORDER_OF_SYSTEM, 0.0 );
      p_Q_INIT( xw[0], xw[1], Qw );
      // compute to the north
      s[ 0 ] = 0.0;
      s[ 1 ] = 1.0;
      DenseVector<double> xn( eb -> get_x( s ) );
      DenseVector<double> Qn( ORDER_OF_SYSTEM, 0.0 );
      p_Q_INIT( xn[0], xn[1], Qn );
      // compute to the north
      s[ 0 ] = 0.0;
      s[ 1 ] = -1.0;
      DenseVector<double> xs( eb -> get_x( s ) );
      DenseVector<double> Qs( ORDER_OF_SYSTEM, 0.0 );
      p_Q_INIT( xs[0], xs[1], Qs );
      // difference for the slopes
      eb -> set_slope_x( ( Qe - Qw ) / ( xe[0] - xw[0] ) );
      eb -> set_slope_y( ( Qn - Qs ) / ( xn[1] - xs[1] ) );
      // set the mid value
      eb -> set_Q_mid( ( Qe + Qw + Qn + Qs ) / 4 );
      ++eb;
    }

#ifdef DEBUG
    std::cout << "DEBUG: mesh constructor complete \n";
#endif

    // now all that remains is to initialise the starting (black) mesh
    eb = BLACK_ELTS.begin();
    while ( eb != BLACK_ELTS.end() )
    {
      DenseVector<double> x( 2, 0.0 );
      DenseVector<double> s( 2, 0.0 );
      x = eb -> get_x( s );
      DenseVector<double> Q( ORDER_OF_SYSTEM, 0.0 );
      p_Q_INIT( x[0], x[1], Q );
      eb -> set_Q_mid( Q );
      ++eb;
    }

    calc_slopes( &BLACK_ELTS );

  }

  TwoD_TVDLF_Mesh::~TwoD_TVDLF_Mesh()
  {}

  void TwoD_TVDLF_Mesh::dump_gnu( std::string filename )
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 6 );
    DenseVector<double> s( 2, 0.0 );
    {
      elt_iter e( BLACK_ELTS.begin() );
      while ( e != BLACK_ELTS.end() )
      {
        dump << e -> get_x( s )[0] << " " << e -> get_x( s )[1] << " ";
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          dump << e -> get_Q( s )[i] << " ";
        }
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          dump << e -> get_slope_x()[i] << " ";
        }
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          dump << e -> get_slope_y()[i] << " ";
        }
        dump << "\n";
        if ( e -> face_is_external( 1 ) )
          dump << "\n";
        ++e;
      }
    }
    dump.close();
  }

  void TwoD_TVDLF_Mesh::dump_nodes_x( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 6 );
    DenseVector<double> s( 2, 0.0 );
    celt_iter e;
    {
      e = BLACK_ELTS.begin();
      for ( unsigned i = 0; i < NX - 1; ++i )
      {
        dump << e -> get_x( s )[0] << "\n";
        ++e;
      }
    }
    dump.close();
  }

  void TwoD_TVDLF_Mesh::dump_nodes_y( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 6 );
    DenseVector<double> s( 2, 0.0 );
    celt_iter e;
    {
      e = BLACK_ELTS.begin();
      for ( unsigned i = 0; i < NY - 1; ++i )
      {
        dump << e -> get_x( s )[1] << "\n";
        e += NX - 1;
      }
    }
    dump.close();
  }

  void TwoD_TVDLF_Mesh::dump_data( std::string filename )
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 6 );
    elt_iter e;
    DenseVector<double> s( 2, 0.0 );
    {
      e = BLACK_ELTS.begin();
      while ( e != BLACK_ELTS.end() )
      {
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          dump << e -> get_Q( s )[i] << " ";
        }
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          dump << e -> get_slope_x()[i] << " ";
        }
        for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
        {
          dump << e -> get_slope_y()[i] << " ";
        }
        dump << "\n";
        ++e;
      }
    }
    dump.close();
  }

  DenseVector<double> TwoD_TVDLF_Mesh::get_point_values( const DenseVector<double>& x )
  {
    elt_iter e( get_elt_iter_from_x( x ) );
    DenseVector<double> s( e-> get_s( x ) );
    return e -> get_Q( s );
  }

  void TwoD_TVDLF_Mesh::set_limiter( const unsigned& id )
  {
    LIMITER = id;
  }

  double TwoD_TVDLF_Mesh::update( const double& CFL, const double& max_dt )
  {
    // integrate the black mesh data onto the red mesh
    {
      elt_iter e = RED_ELTS.begin();
      while ( e != RED_ELTS.end() )
      {
        e -> set_Q_mid( e -> contributed_Q() / e -> get_dA() );
        ++e;
      }
    }
    calc_slopes( &RED_ELTS );

    // determine the first time step from the CFL constraint
    double first_dt;
    {
      elt_iter e = BLACK_ELTS.begin();
      first_dt = e -> get_max_dt();
      ++e;
      while ( e != BLACK_ELTS.end() )
      {
        first_dt = std::min( first_dt, e -> get_max_dt() );
        ++e;
      }
      first_dt *= CFL;
    }
    if ( first_dt > max_dt / 2 )
    {
      first_dt = max_dt / 2;
    }

    actions_before_time_step1( first_dt );
    // do the time step
    {
      elt_iter e = RED_ELTS.begin();
      while ( e != RED_ELTS.end() )
      {
        // add to the corrections
        e -> add_flux_contributions( first_dt );
        ++e;
      }
    }
    calc_slopes( &RED_ELTS );
    MESH_TIME += first_dt;

    // integrate the red mesh data onto the black mesh
    {
      elt_iter e = BLACK_ELTS.begin();
      while ( e != BLACK_ELTS.end() )
      {
        e -> set_Q_mid( e -> contributed_Q() / e -> get_dA() );
        ++e;
      }
    }
    calc_slopes( &BLACK_ELTS );

    // determine the first time step from the CFL constraint
    double second_dt;
    {
      elt_iter e = RED_ELTS.begin();
      second_dt = e -> get_max_dt();
      ++e;
      while ( e != RED_ELTS.end() )
      {
        second_dt = std::min( second_dt, e -> get_max_dt() );
        ++e;
      }
      second_dt *= CFL;
    }
    if ( first_dt + second_dt > max_dt )
    {
      second_dt = max_dt - first_dt;
    }

    actions_before_time_step2( second_dt );
    // do the time step
    {
      elt_iter e = BLACK_ELTS.begin();
      while ( e != BLACK_ELTS.end() )
      {
        // add to the corrections
        e -> add_flux_contributions( second_dt );
        ++e;
      }
    }
    calc_slopes( &BLACK_ELTS );
    MESH_TIME += second_dt;

    //std::cout << "    currently at " << MESH_TIME << "\n";

    return first_dt + second_dt;
  }

  void TwoD_TVDLF_Mesh::update_to( const double& CFL, const double& t_end )
  {
    do
    {
      std::cout << "    computing to " << t_end << "\n";
      update( CFL, std::abs( t_end - MESH_TIME ) );
    }
    while ( MESH_TIME < t_end );
  }

  const double& TwoD_TVDLF_Mesh::get_time() const
  {
    return MESH_TIME;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::integrate( std::string mesh_colour )
  {
    vector_of_elts* elts( get_elts_from_colour( mesh_colour ) );
    elt_iter e = elts -> begin();
    DenseVector<double> I( ORDER_OF_SYSTEM, 0.0 );
    while ( e != elts -> end() )
    {
      const DenseVector<double> sw( 2, -1.0 );
      const DenseVector<double> ne( 2, 1.0 );
      I += e -> get_int_Q( sw, ne );
      ++e;
    }
    return I;
  }

  std::vector<TwoD_TVDLF_Elt>::iterator TwoD_TVDLF_Mesh::get_elt_iter_from_x( const DenseVector<double>& x, std::string mesh_colour )
  {
    elt_iter found_elt;
    vector_of_elts* elts( get_elts_from_colour( mesh_colour ) );
    std::size_t Mx( get_number_elts_in_x( mesh_colour ) );

    elt_iter e = elts -> begin();
    while ( e != elts -> end() )
    {
      // local coordinates of the NE and SW corners
      DenseVector<double> ne( 2, 1.0 );
      DenseVector<double> sw( -ne );
      double west = e -> get_x( sw )[ 0 ];
      double east = e -> get_x( ne )[ 0 ];
      if ( ( x[ 0 ] >= west ) && ( x[ 0 ] <= east ) )
      {
        break;
      }
      ++e;
    }

    while ( true )
    {
      // local coordinates of the NE and SW corners
      DenseVector<double> ne( 2, 1.0 );
      DenseVector<double> sw( -ne );
      double south = e -> get_x( sw )[ 1 ];
      double north = e -> get_x( ne )[ 1 ];
      if ( ( x[ 1 ] >= south ) && ( x[ 1 ] <= north ) )
      {
        break;
      }
      if ( e + Mx < elts -> end() )
      {
        e += Mx;
      }
      else
      {
        // if we are here, then we're looking for a point outside the mesh
        std::string problem;
        problem = "The TwoD_TVDLF_Mesh::get_elt_iter_from_x method has been called\n";
        problem += "for a position that is outside the boundaries of the mesh.\n";
        throw ExceptionRuntime( problem );
      }
    }
    // if we get to here, then we've found the elt.
    return e;
  }


  std::vector<TwoD_TVDLF_Elt>* TwoD_TVDLF_Mesh::get_elts_from_colour( std::string mesh_colour )
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

  std::size_t TwoD_TVDLF_Mesh::get_number_elts_in_x( std::string mesh_colour )
  {
    std::size_t N;
    if ( mesh_colour == "black" )
    {
      N = NX - 1;
    }
    else
      if ( mesh_colour == "red" )
      {
        N = NX;
      }
      else
      {
        std::string problem;
        problem = " The TwoD_TVDLF_Mesh object has been passed an unrecognised mesh identifier. \n";
        problem += " valid options are 'black' or 'red'.\n";
        throw ExceptionRuntime( problem );
      }
    return N;
  }

  void TwoD_TVDLF_Mesh::calc_slopes( vector_of_elts* elt_vector )
  {
    set_boundary_Q( elt_vector );
    elt_iter e( elt_vector -> begin() );
    DenseVector<double> slope_east( ORDER_OF_SYSTEM, 0.0 );
    DenseVector<double> slope_west( ORDER_OF_SYSTEM, 0.0 );
    DenseVector<double> slope_north( ORDER_OF_SYSTEM, 0.0 );
    DenseVector<double> slope_south( ORDER_OF_SYSTEM, 0.0 );
    switch ( LIMITER )
    {
    case 0:
      // minmod
      while ( e != elt_vector -> end() )
      {
        slope_east = east_diff( e );
        slope_west = west_diff( e );
        slope_south = south_diff( e );
        slope_north = north_diff( e );
        e -> set_slope_x( minmod( slope_east, slope_west ) );
        e -> set_slope_y( minmod( slope_north, slope_south ) );
        ++e;
      }
      break;
    case 1:
      // MC
      while ( e != elt_vector -> end() )
      {
        slope_east = east_diff( e );
        slope_west = west_diff( e );
        DenseVector<double> slope_ew( EW_diff( e ) );
        slope_south = south_diff( e );
        slope_north = north_diff( e );
        DenseVector<double> slope_ns( NS_diff( e ) );
        const DenseVector<double> temp_x( minmod( slope_east * 2, slope_west * 2 ) );
        const DenseVector<double> slope_x( minmod( temp_x, slope_ew ) );
        const DenseVector<double> temp_y( minmod( slope_north * 2, slope_south * 2 ) );
        const DenseVector<double> slope_y( minmod( temp_y, slope_ns ) );
        e -> set_slope_x( slope_x );
        e -> set_slope_y( slope_y );
        ++e;
      }
    case 2:
      /// superbee
      while ( e != elt_vector -> end() )
      {
        slope_east = east_diff( e );
        slope_west = west_diff( e );
        slope_south = south_diff( e );
        slope_north = north_diff( e );
        const DenseVector<double> slope_x = maxmod(
                                              minmod( slope_east, slope_west * 2. ),
                                              minmod( slope_east * 2., slope_west )
                                            );
        e -> set_slope_x( slope_x );
        const DenseVector<double> slope_y = maxmod(
                                              minmod( slope_north, slope_south * 2. ),
                                              minmod( slope_north * 2., slope_south )
                                            );
        e -> set_slope_y( slope_y );
        ++e;
      }
      break;
    default:
      std::string problem;
      problem = " The TwoD_TVDLF_Mesh object has an unrecognised 'LIMITER' identifier. \n";
      throw ExceptionRuntime( problem );
    }
  }

  DenseVector<double> TwoD_TVDLF_Mesh::east_diff( elt_iter e ) const
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    std::set< int > faces( e -> get_external_faces() );
    const int face_index( 1 );
    if ( e -> face_is_external( face_index ) )
    {
      // interior difference
      diff = west_diff( e );
      boundary_diff( e, face_index, diff );
    }
    else
    {
      // We're not on an edge, so we can difference
      elt_iter j( e -> get_ptrs( face_index ) );
      diff = ( j -> get_Q_mid() - e -> get_Q_mid() )
             / ( j -> get_x_mid()[0] - e -> get_x_mid()[0] );
    }
    return diff;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::west_diff( elt_iter e ) const
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    std::set< int > faces( e -> get_external_faces() );
    const int face_index( 3 );
    if ( e -> face_is_external( face_index ) )
    {
      // interior difference
      diff = east_diff( e );
      boundary_diff( e, face_index, diff );
    }
    else
    {
      // We're not on an edge, so we can difference
      elt_iter j( e -> get_ptrs( face_index ) );
      diff = ( e -> get_Q_mid() - j -> get_Q_mid() )
             / ( e -> get_x_mid()[0] - j -> get_x_mid()[0] );
    }
    return diff;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::north_diff( elt_iter e ) const
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    std::set< int > faces( e -> get_external_faces() );
    const int face_index( 2 );
    if ( e -> face_is_external( face_index ) )
    {
      // interior difference
      diff = south_diff( e );
      boundary_diff( e, face_index, diff );
    }
    else
    {
      // We're not on an edge, so we can difference
      elt_iter j( e -> get_ptrs( face_index ) );
      diff = ( j -> get_Q_mid() - e -> get_Q_mid() )
             / ( j -> get_x_mid()[1] - e -> get_x_mid()[1] );
    }
    return diff;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::south_diff( elt_iter e ) const
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    std::set< int > faces( e -> get_external_faces() );
    const int face_index( 0 );
    if ( e -> face_is_external( face_index ) )
    {
      // interior difference
      diff = north_diff( e );
      boundary_diff( e, face_index, diff );
    }
    else
    {
      // We're not on an edge, so we can difference
      elt_iter j( e -> get_ptrs( face_index ) );
      diff = ( e -> get_Q_mid() - j -> get_Q_mid() )
             / ( e -> get_x_mid()[1] - j -> get_x_mid()[1] );
    }
    return diff;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::NS_diff( elt_iter e )
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    std::set< int > faces( e -> get_external_faces() );
    if ( e -> face_is_external( 0 ) )
    {
      // interior difference
      diff = north_diff( e );
      boundary_diff( e, 0, diff );
    }
    else
      if ( e -> face_is_external( 2 ) )
      {
        // interior difference
        diff = south_diff( e );
        boundary_diff( e, 2, diff );
      }
      else
      {
        // We're not on an edge, so we can difference
        elt_iter jS( e -> get_ptrs( 0 ) );
        elt_iter jN( e -> get_ptrs( 2 ) );
        diff = ( jN -> get_Q_mid() - jS -> get_Q_mid() )
               / ( jN -> get_x_mid()[1] - jS -> get_x_mid()[1] );
      }
    return diff;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::EW_diff( elt_iter e )
  {
    DenseVector<double> diff( ORDER_OF_SYSTEM, 0.0 );
    std::set< int > faces( e -> get_external_faces() );
    if ( e -> face_is_external( 1 ) )
    {
      // interior difference
      diff = west_diff( e );
      boundary_diff( e, 1, diff );
    }
    else
      if ( e -> face_is_external( 3 ) )
      {
        // interior difference
        diff = east_diff( e );
        boundary_diff( e, 3, diff );
      }
      else
      {
        // We're not on an edge, so we can difference
        elt_iter jE( e -> get_ptrs( 1 ) );
        elt_iter jW( e -> get_ptrs( 3 ) );
        diff = ( jE -> get_Q_mid() - jW -> get_Q_mid() )
               / ( jE -> get_x_mid()[0] - jW -> get_x_mid()[0] );
      }
    return diff;
  }

  int TwoD_TVDLF_Mesh::sgn( double a ) const
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

  DenseVector<double> TwoD_TVDLF_Mesh::minmod( DenseVector<double> A, DenseVector<double> B ) const
  {
    DenseVector<double> MM;
    for ( std::size_t i = 0; i < A.size(); ++i )
    {
      MM.push_back( 0.5 * ( sgn( A[i] ) + sgn( B[i] ) )
                    * std::min( std::abs( A[i] ), std::abs( B[i] ) ) );
    }
    return MM;
  }

  DenseVector<double> TwoD_TVDLF_Mesh::maxmod( DenseVector<double> A, DenseVector<double> B ) const
  {
    DenseVector<double> MM;
    for ( std::size_t i = 0; i < A.size(); ++i )
    {
      MM.push_back( 0.5 * ( sgn( A[i] ) + sgn( B[i] ) )
                    * std::max( std::abs( A[i] ), std::abs( B[i] ) ) );
    }
    return MM;
  }

  void TwoD_TVDLF_Mesh::boundary_diff( elt_iter e, const int& face_index, DenseVector<double>& diff ) const
  {
    DenseVector<double> se( 2, 0.0 );
    // get the appropriate local coordinate for the mid point on the edge
    switch ( face_index )
    {
    case 0:
      // south
      se[ 0 ] = 0.0;
      se[ 1 ] = -1.0;
      break;
    case 1:
      // east
      se[ 0 ] = 1.0;
      se[ 1 ] = 0.0;
      break;
    case 2:
      // north
      se[ 0 ] = 0.0;
      se[ 1 ] = 1.0;
      break;
    case 3:
      // west
      se[ 0 ] = -1.0;
      se[ 1 ] = 0.0;
      break;
    }
    // for the edge value
    DenseVector<double> Q( e -> get_Q( se ) );
    // look at the user-defined BCs for any defined Dirichlet conditions
    std::vector<bool> inflow = e -> p_system -> edge_values( face_index, e -> get_x( se ), Q );
    // the default is a zero slope in accordance with Levy & Tadmor (1997)
    DenseVector<double> sigma_n( ORDER_OF_SYSTEM, 0.0 );
    // allow the user to override the zero slope
    e -> p_system -> edge_slopes( face_index, e -> get_x( se ), sigma_n );
    // use edge values for inflow conditions
    for ( std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i )
    {
      if ( inflow[ i ] == true )
      {
        // set a zero slope as per Levy & Tadmor (1997)?
        diff[ i ] = sigma_n[ i ];// ( Qe[ i ] - Qm[ i ] ) / delta;
        //std::cout << face_index << " " << diff[ i ] << "\n";
      }
    }
  }

}
