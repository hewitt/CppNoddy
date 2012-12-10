/// \file TwoD_Node_Mesh.cpp
/// Implementation of a two dimensional mesh object. Data
/// is stored on a nodal mesh.

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <Exceptions.h>
#include <DenseVector.h>
#include <DenseMatrix.h>
#include <TwoD_Node_Mesh.h>
#include <Utility.h>

namespace CppNoddy
{

  template <typename _Type>
  void TwoD_Node_Mesh<_Type>::set_nodes_vars( const std::size_t nodex, const std::size_t nodey, const DenseVector<_Type>& U )
  {
#ifdef PARANOID
    if ( U.size() > NV )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.set_nodes_vars method is trying to use a \n";
      problem += " vector that has more entries than variables stored in the mesh. \n";
      throw ExceptionRuntime( problem );
    }
#endif
    // assign contents of U to the member data
    std::size_t offset( ( nodex * NY + nodey ) * NV );
    for ( std::size_t var = 0; var < NV; ++var )
    {
      VARS[ offset++ ] = U[ var ];
    }
  }

  template <typename _Type>
  DenseVector<_Type> TwoD_Node_Mesh<_Type>::get_nodes_vars( const std::size_t nodex, const std::size_t nodey ) const
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.get_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
#endif
    // construct a vector with NV elements starting from a pointer
    DenseVector<_Type> nodes_vars( NV, &VARS[ ( nodex * NY + nodey ) * NV ] );
    return nodes_vars;
  }

  template <typename _Type>
  void TwoD_Node_Mesh<_Type>::assign( const _Type elt )
  {
    VARS.assign( NX * NY * NV, elt );
  }

  template <typename _Type>
  std::pair< std::size_t, std::size_t > TwoD_Node_Mesh<_Type>::get_nnodes() const
  {
    std::pair< std::size_t, std::size_t > nodes;
    nodes.first = NX;
    nodes.second = NY;
    return nodes;
  }

  template <typename _Type>
  std::size_t TwoD_Node_Mesh<_Type>::get_nvars() const
  {
    return NV;
  }

  template <typename _Type>
  const DenseVector<double>& TwoD_Node_Mesh<_Type>::xnodes() const
  {
    return X;
  }

  template <typename _Type>
  const DenseVector<double>& TwoD_Node_Mesh<_Type>::ynodes() const
  {
    return Y;
  }

  template <typename _Type>
  DenseMatrix<_Type> TwoD_Node_Mesh<_Type>::get_var_as_matrix( std::size_t var ) const
  {
#ifdef PARANOID
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.get_var_as_matrix method is trying to use a \n";
      problem += " variable index bigger than the number of variables in the mesh. \n";
      throw ExceptionRange( problem, NV, var );
    }
#endif
    DenseMatrix<_Type> temp( NX, NY, 0.0 );
    for ( std::size_t i = 0; i < NX; ++i )
    {
      for  ( std::size_t j = 0; j < NY; ++j )
      {
        temp( i, j ) = VARS[ ( i * NY + j ) * NV + var ];
      }
    }
    return temp;
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::remesh1( const DenseVector<double>& newX, const DenseVector<double>& newY )
  {
#ifdef PARANOID
    // check start & end 
    if ( std::abs( X[ 0 ] - newX[ 0 ] ) > 1.e-10 ||
         std::abs( X[ X.size() - 1 ] - newX[ newX.size() - 1 ] ) > 1.e-10 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.remesh1 method has been called with \n";
      problem += " a passed X coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw ExceptionRuntime( problem );
    }
    // check monotonic node positions
    for ( std::size_t i = 0; i < newX.size() - 1; ++i )
    {
      if ( newX[ i ] >= newX[ i + 1 ] )
      {
        std::string problem;
        problem = " The TwoD_Node_Mesh.remesh1 method has been passed \n";
        problem += " a non-monotonic X coordinate vector. \n";
        problem += Utility::stringify( newX[ i ], 6 ) + " vs. " + Utility::stringify( newX[ i + 1 ], 6 );
        throw ExceptionRuntime( problem );
      }
    }
    // check start and end
    if ( std::abs( Y[ 0 ] - newY[ 0 ] ) > 1.e-10 ||
         std::abs( Y[ Y.size() - 1 ] - newY[ newY.size() - 1 ] ) > 1.e-10 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.remesh1 method has been called with \n";
      problem += " a passed Y coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw ExceptionRuntime( problem );
    }
    // check monotonic node positions
    for ( std::size_t i = 0; i < newY.size() - 1; ++i )
    {
      if ( newY[ i ] >= newY[ i + 1 ] )
      {
        std::string problem;
        problem = " The TwoD_Node_Mesh.remesh1 method has been passed \n";
        problem += " a non-monotonic Y coordinate vector. \n";
        problem += Utility::stringify( newY[ i ], 6 ) + " vs. " + Utility::stringify( newY[ i + 1 ], 6 );
        throw ExceptionRuntime( problem );
      }
    }
#endif

    // new variables storage
    DenseVector<_Type> newvars( newX.size() * newY.size() * NV, 0.0 );

    // left boundary
    {
      std::size_t xnode( 0 );
      // bottom left corner copy
      for ( unsigned var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + 0 ) * NV + var ] = get_nodes_vars( 0, 0 )[ var ];
      }
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaY( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y.size() - 1; ++j )
        {
          if ( ( Y[ j ] <= newY[ ynode ] ) && ( newY[ ynode ] < Y[ j + 1 ] ) )
          {
            below_j = j;
            deltaY = newY[ ynode ] - Y[ j ];
          }
        }
        DenseVector<_Type> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 ) - get_nodes_vars( left_i, below_j ) )
                                     / ( coord( left_i, below_j + 1 ).second - coord( left_i, below_j ).second );
        DenseVector<_Type> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdY * deltaY;
        for ( unsigned var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
      // top left corner copy
      for ( unsigned var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + newY.size() - 1 ) * NV + var ] = get_nodes_vars( 0, NY - 1 )[ var ];
      }
    }
    // right boundary
    {
      std::size_t xnode( newX.size() - 1 );
      // bottom right corner copy
      for ( unsigned var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + 0 ) * NV + var ] = get_nodes_vars( NX - 1, 0 )[ var ];
      }
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( X.size() - 1 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaY( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y.size() - 1; ++j )
        {
          if ( ( Y[ j ] <= newY[ ynode ] ) && ( newY[ ynode ] < Y[ j + 1 ] ) )
          {
            below_j = j;
            deltaY = newY[ ynode ] - Y[ j ];
          }
        }
        DenseVector<_Type> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 ) - get_nodes_vars( left_i, below_j ) )
                                     / ( coord( left_i, below_j + 1 ).second - coord( left_i, below_j ).second );
        DenseVector<_Type> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdY * deltaY;
        for ( unsigned var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
      // bottom right corner copy
      for ( unsigned var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + newY.size() - 1 ) * NV + var ] = get_nodes_vars( NX - 1, NY - 1 )[ var ];
      }
    }
    // bottom boundary
    {
      std::size_t ynode( 0 );
      for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaX( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X.size() - 1; ++i )
        {
          if ( ( X[ i ] <= newX[ xnode ] ) && ( newX[ xnode ] < X[ i + 1 ] ) )
          {
            left_i = i;
            deltaX = newX[ xnode ] - X[ i ];
          }
        }
        DenseVector<_Type> dvarsdX = ( get_nodes_vars( left_i + 1, below_j ) - get_nodes_vars( left_i, below_j ) )
                                     / ( coord( left_i + 1, below_j ).first - coord( left_i, below_j ).first );
        DenseVector<_Type> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdX * deltaX;
        for ( unsigned var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
    }
    // top boundary
    {
      std::size_t ynode( newY.size() - 1 );
      for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( Y.size() - 1 ); // bracketing index
        double deltaX( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X.size() - 1; ++i )
        {
          if ( ( X[ i ] <= newX[ xnode ] ) && ( newX[ xnode ] < X[ i + 1 ] ) )
          {
            left_i = i;
            deltaX = newX[ xnode ] - X[ i ];
          }
        }
        DenseVector<_Type> dvarsdX = ( get_nodes_vars( left_i + 1, below_j ) - get_nodes_vars( left_i, below_j ) )
                                     / ( coord( left_i + 1, below_j ).first - coord( left_i, below_j ).first );
        DenseVector<_Type> interpolated_vars = get_nodes_vars( left_i, below_j ) + dvarsdX * deltaX;
        for ( unsigned var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
    }
    // loop thru interior nodes of the destination mesh one node at a time
    for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
    {
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X.size() - 1; ++i )
        {
          if ( ( X[ i ] <= newX[ xnode ] ) && ( newX[ xnode ] < X[ i + 1 ] ) )
          {
            left_i = i;
          }
        }
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y.size() - 1; ++j )
        {
          if ( ( Y[ j ] <= newY[ ynode ] ) && ( newY[ ynode ] < Y[ j + 1 ] ) )
          {
            below_j = j;
          }
        }
        DenseVector<_Type> dvarsdX = ( get_nodes_vars( left_i + 1, below_j ) - get_nodes_vars( left_i, below_j ) )
                                     / ( coord( left_i + 1, below_j ).first - coord( left_i, below_j ).first );
        DenseVector<_Type> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 ) - get_nodes_vars( left_i, below_j ) )
                                     / ( coord( left_i, below_j + 1 ).second - coord( left_i, below_j ).second );

        DenseVector<_Type> interpolated_vars_bottom =
          ( get_nodes_vars( left_i, below_j ) * ( coord( left_i + 1, below_j ).first - newX[ xnode ] )
            + get_nodes_vars( left_i + 1, below_j ) * ( newX[ xnode ] - coord( left_i, below_j ).first ) ) /
          ( coord( left_i + 1, below_j ).first - coord( left_i, below_j ).first );

        DenseVector<_Type> interpolated_vars_top =
          ( get_nodes_vars( left_i, below_j + 1 ) * ( coord( left_i + 1, below_j + 1 ).first - newX[ xnode ] )
            + get_nodes_vars( left_i + 1, below_j + 1 ) * ( newX[ xnode ] - coord( left_i, below_j + 1 ).first ) ) /
          ( coord( left_i + 1, below_j + 1 ).first - coord( left_i, below_j + 1 ).first );

        DenseVector<_Type> interpolated_vars =
          (  interpolated_vars_bottom * ( coord( left_i, below_j + 1 ).second - newY[ ynode ] )
             +  interpolated_vars_top * ( newY[ ynode ] - coord( left_i, below_j ).second ) ) /
          ( coord( left_i, below_j + 1 ).second - coord( left_i, below_j ).second );

        for ( unsigned var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] = interpolated_vars[ var ];
        }
      }
    }
    // finally replace the old nodes with the new ones
    X = newX;
    Y = newY;
    NX = newX.size();
    NY = newY.size();
    VARS = newvars;
  }

  template<typename _Type>
  OneD_Node_Mesh<_Type> TwoD_Node_Mesh<_Type>::get_xsection_at_xnode( const std::size_t nodex ) const
  {
    OneD_Node_Mesh<_Type> xsection( Y, NV );
    for ( std::size_t nodey = 0; nodey < NY; ++nodey )
    {
      xsection.set_nodes_vars( nodey, this -> get_nodes_vars( nodex, nodey ) );
    }
    return xsection;
  }

  template<typename _Type>
  OneD_Node_Mesh<_Type> TwoD_Node_Mesh<_Type>::get_xsection_at_ynode( const std::size_t nodey ) const
  {
    OneD_Node_Mesh<_Type> xsection( X, NV );
    for ( std::size_t nodex = 0; nodex < NX; ++nodex )
    {
      xsection.set_nodes_vars( nodex, this -> get_nodes_vars( nodex, nodey ) );
    }
    return xsection;
  }


  template <typename _Type>
  void TwoD_Node_Mesh<_Type>::dump() const
  {
    for ( std::size_t var = 0; var < NV; ++var )
    {
      std::cout << "Variable : " << var << "\n";
      std::cout << " x = ";
      for ( std::size_t i = 0; i < NX; ++i )
      {
        std::cout << X[ i ] << ", ";
      }
      std::cout << "\n";
      for ( std::size_t j = 0; j < NY; ++j )
      {
        std::cout << " y = " << Y[ j ] << "\n";
        for ( std::size_t i = 0; i < NX; ++i )
        {
          std::cout << VARS[ ( i * NY + j ) * NV + var ] << ", ";
        }
        std::cout << "\n";
      }
    }
  }

  template<>
  void TwoD_Node_Mesh<double>::dump_gnu( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );

    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        dump << X[ i ] << " " << Y[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

 template <>
 void TwoD_Node_Mesh<D_complex>::dump_gnu( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );

    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        dump << X[ i ] << " " << Y[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << real(VARS[ ( i * NY + j ) * NV + var ]) << " ";
          dump << imag(VARS[ ( i * NY + j ) * NV + var ]) << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::dump_var( std::string filename, const unsigned var ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    dump.precision( 9 );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        dump << VARS[ ( i * NY + j ) * NV + var ] << "\n";
      }
    }
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::dump( std::string filename ) const
  {
    std::ofstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    //dump << NX << " " << NY << " " << NV << "\n";
    dump.precision( 9 );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        dump << X[ i ] << " " << Y[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        dump << "\n";
      }
    }
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::read( std::string filename, bool reset )
  {
    std::ifstream dump;
    dump.open( filename.c_str() );
    dump.precision( 9 );
    dump.setf( std::ios::showpoint );
    dump.setf( std::ios::showpos );
    dump.setf( std::ios::scientific );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        double x, y;
        dump >> x;
        dump >> y;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          double value;
          dump >> value;
          VARS[ ( i * NY + j ) * NV + var ] = value;
        }
        if ( reset != true )
        {
          // if not reseting the mesh we should check the node positions
          if ( ( std::abs( x - X[ i ] ) > 1.e-6 ) || ( std::abs( y - Y[ j ] ) > 1.e-6 ) )
          {
            std::cout << " Read x = " << x << " Expected x = " << X[ i ] << "; Read y = " << y << " Expected y = " << Y[ j ] << " \n";
            std::cout << " Absolute differences are " << abs( x - X[i] ) << " and " << abs( y - Y[j] ) << "\n";              
            std::string problem;
            problem = " The TwoD_Node_Mesh.read method is trying to read a \n";
            problem += " file whose nodal points are in a different position. \n";
            throw ExceptionRuntime( problem );
          }
        }
        else
        {
          X[ i ] = x;
          Y[ j ] = y;
        }
      }
    }
  }

  //the templated versions we require are:
  template class TwoD_Node_Mesh<double>
  ;
  template class TwoD_Node_Mesh<std::complex<double> >
  ;

}
