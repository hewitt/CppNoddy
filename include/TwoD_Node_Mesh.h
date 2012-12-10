/// \file TwoD_Node_Mesh.h
/// A specification for a two dimensional mesh object. Data
/// is stored on a nodal mesh.
#ifndef TWOD_NODE_MESH_H
#define TWOD_NODE_MESH_H

#include <vector>
#include <fstream>

#include <DenseVector.h>
#include <DenseMatrix.h>
#include <OneD_Node_Mesh.h>

namespace CppNoddy
{

  /// A two dimensional mesh utility object.
  template <typename _Type>
  class TwoD_Node_Mesh
  {
  public:

    TwoD_Node_Mesh()
    {}

    /// ctor
    TwoD_Node_Mesh( const DenseVector<double>& x_nodes, const DenseVector<double>& y_nodes, const std::size_t nvars ) :
        NX( x_nodes.size() ), NY( y_nodes.size() ), NV( nvars ), X( x_nodes ), Y( y_nodes )
    {
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      VARS = DenseVector<_Type>( NX * NY * NV, 0.0 );
    }

    // ctor from a file
    TwoD_Node_Mesh( std::string filename ) 
    {
      read( filename, true );   
    }

    /// dtor
    virtual ~TwoD_Node_Mesh()
    {}

    /// Access operator for a nodal point/variable in the mesh
    /// \param nodex The nodal index value in the first direction
    /// \param nodey The nodal index value in the second direction
    /// \param var The variable index to be accessed
    _Type& operator()( const std::size_t nodex, const std::size_t nodey, const std::size_t var );

    /// Const access operator for a nodal point/variable in the mesh
    /// \param nodex The nodal index value in the first direction
    /// \param nodey The nodal index value in the second direction
    /// \param var The variable index to be accessed
    const _Type& operator()( const std::size_t nodex, const std::size_t nodey, const std::size_t var ) const;

    /// Access the nodal position - to ensure the uniformity
    /// of the mesh nodal positions, we do not allow any
    /// non-const access to the X, Y vectors in the base class.
    /// \param nodex The x nodal position to return
    /// \param nodey The y nodal position to return
    /// \return The spatial position of this node as a pair
    std::pair<double, double> coord( const std::size_t nodex, const std::size_t nodey ) const;

    /// Set the variables stored at A SPECIFIED node
    /// \param nodex The x nodal index to be set
    /// \param nodey The y nodal index to be set
    /// \param U The vector of VARIABLES to be written to this nodal point
    void set_nodes_vars( const std::size_t nodex, const std::size_t nodey, const DenseVector<_Type>& U );

    /// Get the variables stored at A SPECIFIED node
    /// \param nodex The x nodal index to be returned
    /// \param nodey The y nodal index to be returned
    /// \return The vector of VARIABLES stored at this nodal point
    DenseVector<_Type> get_nodes_vars( const std::size_t nodex, const std::size_t nodey ) const;

    /// Get a cross section of the 2D mesh at a specified (constant) x node
    /// \param nodex The x nodal index at which the cross section is to be taken
    /// \return A 1D nodal mesh  
    OneD_Node_Mesh<_Type> get_xsection_at_xnode( const std::size_t nodex ) const;

    /// Get a cross section of the 2D mesh at a specified (constant) y node
    /// \param nodey The y nodal index at which the cross section is to be taken
    /// \return A 1D nodal mesh  
    OneD_Node_Mesh<_Type> get_xsection_at_ynode( const std::size_t nodey ) const; 

    /// Assign an element to all entries in the mesh
    /// \param elt The element to be assigned to the mesh
    void assign( const _Type elt );

    /// Get the number of nodes in the two directions of the 2D mesh
    /// \return A pair consisting of the number of nodes in the 2 directions
    std::pair< std::size_t, std::size_t > get_nnodes() const;

    /// Get the number of variables that are stored at each node
    /// \return The number of variables that have data stored at
    /// each nodal point
    std::size_t get_nvars() const;

    /// Access the vector of x-nodal positions
    /// \return A vector of the nodal positions for this mesh
    const DenseVector<double>& xnodes() const;

    /// Access the vector of y-nodal positions
    /// \return A vector of the nodal positions for this mesh
    const DenseVector<double>& ynodes() const;

    /// Return a matrix corresponding to each nodal point in the mesh
    /// Each matrix element will contain a specified variable number
    /// \param var The variable number to be accessed
    /// \return A dense matrix of the specified variable
    DenseMatrix<_Type> get_var_as_matrix( std::size_t var ) const;

    /// Interpolate this mesh data (bilinearly) into a new
    /// mesh with nodal points defined in the argument list.
    /// Not written to be efficient, so you probably don't want
    /// to do any repeated calls with this method.
    /// \param newX The x-nodal coordinates to be used in the new mesh.
    /// \param newY The y-nodal coordinates to be used in the new mesh.
    void remesh1( const DenseVector<double>& newX, const DenseVector<double>& newY );

    /// A simple method for dumping data to std::cout
    void dump() const;

    /// A simple method for dumping data to a file
    /// \param filename The filename to write the data to (will overwrite)
    void dump( std::string filename ) const;
    
    /// A simple method for dumping a single variable to a file with no nodal information
    /// \param filename The filename to write the data to (will overwrite)
    /// \param var The index of the variable to be dumped to output
    void dump_var( std::string filename, const unsigned var ) const;

    /// A simple method for reading data from a file
    /// \param filename The filename to write the data to (will overwrite)
    /// \param reset Will reset the nodal positions using those from the file  
    void read( std::string filename, const bool reset = false );

    /// A simple method for dumping data to a file for gnuplot surface plotting
    /// \param filename The filename to write the data to (will overwrite)
    void dump_gnu( std::string filename ) const;


    /// Get a bilinearly interpolated value at a specified point
    /// \param x x-coordinate in the 2D mesh      
    /// \param y y-coordinate in the 2D mesh      
    /// \return A vector of bilinearly interpolated values
    DenseVector<_Type> get_interpolated_vars( const double& x, const double& y)
    {
      // check start and end
      if ( ( x < X[0] ) || ( x>X[NX-1] ) )
      {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method has been called with \n";
        problem += " an x coordinate that lies outside the mesh. \n";
        throw ExceptionRuntime( problem );
      }
      // check start and end
      if ( ( y < Y[0] ) || ( y>Y[NY-1] ) )
      {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method has been called with \n";
        problem += " a y coordinate that lies outside the mesh. \n";
        throw ExceptionRuntime( problem );
      }
      int bottom_j(-2);
      for ( unsigned j = 0; j < NY-1; ++j )
	{
	  if ( ( y > Y[j] ) && ( y < Y[j+1] ) )
	    {
	      bottom_j = j;
	    }
	  if ( ( abs(y-Y[j])<1.e-10 ) || ( abs(y-Y[j+1])<1.e-10 ) )
	    {
	      bottom_j = j;
	    }
	}
      if ( bottom_j == -1 )
	{
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method is broken.\n";
        throw ExceptionRuntime( problem );
	}

      std::cout << y << " " << Y[bottom_j] << " " << Y[bottom_j+1] << "\n";
      //
      OneD_Node_Mesh<_Type> bottom_row = get_xsection_at_ynode( bottom_j );
      OneD_Node_Mesh<_Type> top_row = get_xsection_at_ynode( bottom_j+1 );
      const double y1 = Y[ bottom_j ]; 
      const double y2 = Y[ bottom_j+1 ];
      DenseVector<_Type> result = top_row.get_interpolated_vars(x)*( y2-y )/( y2-y1 )
        + bottom_row.get_interpolated_vars(x)*( y-y1 )/( y2-y1 );
      std::cout << "x,y,interp: " << x << " " << y << " " << result[0] << "\n"; 
      return result; 
    }
    
            
  protected:

    // we'll store the number of nodes
    std::size_t NX, NY, NV;
    // store x nodal points
    DenseVector<double> X;
    // store y nodal points
    DenseVector<double> Y;
    // store the nodal values
    ///std::vector< DenseMatrix<_Type> > VARS;
    DenseVector<_Type> VARS;
  };

  template <typename _Type>
  inline _Type& TwoD_Node_Mesh<_Type>::operator()( const std::size_t nodex, const std::size_t nodey, const std::size_t var )
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange( problem, NV, var );
    }
#endif
    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <typename _Type>
  inline const _Type& TwoD_Node_Mesh<_Type>::operator()( const std::size_t nodex, const std::size_t nodey, const std::size_t var ) const
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange( problem, NV, var );
    }
#endif
    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <typename _Type>
  inline std::pair<double, double> TwoD_Node_Mesh<_Type>::coord( const std::size_t nodex, const std::size_t nodey ) const
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Node_Mesh.coord method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
#endif
    std::pair< double, double > pos;
    pos.first = X[ nodex ];
    pos.second = Y[ nodey ];
    return pos;
  }

}

#endif // TWOD_NODE_MESH_H
