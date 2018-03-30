/// \file TwoD_Mapped_Node_Mesh.h
/// A specification for a two dimensional mapped mesh object. Data
/// is stored on a nodal mesh that is non-uniform (in general). This
/// class essentially combines the mapping function for the mesh into
/// the 2D storage of the data.

#ifndef TWOD_MAPPED_NODE_MESH_H
#define TWOD_MAPPED_NODE_MESH_H

#include <vector>
#include <fstream>

#include <DenseVector.h>
#include <Types.h>

namespace CppNoddy
{

  /// A two dimensional (mapped) mesh utility object.
  template <typename _Type>
    class TwoD_Mapped_Node_Mesh
  {
  public:

    /// Mapping function that provides a computational X coordinate
    /// from a physical coordinate.
    /// \param zeta The physical coordinate
    /// \return The corresponding computational coordinate
    virtual double FnComp_X( const double& zeta ) const
    {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_X method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime( problem );
    }
    
    /// Mapping function that provides the first derivative of the
    /// computational X coordinate as a function of the physical coordinate.
    /// \param zeta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Xd( const double& zeta ) const
    {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Xd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime( problem );
    }

    /// Mapping function that provides the second derivative of the
    /// computational X coordinate as a function of the physical coordinate.
    /// \param zeta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Xdd( const double& zeta ) const
    {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Xdd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime( problem );
    }

    /// Mapping function that provides the 
    /// computational Y coordinate as a function of the physical coordinate.
    /// \param eta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Y( const double& eta ) const
    {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Y method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime( problem );
    }

    /// Mapping function that provides the first derivative of the
    /// computational Y coordinate as a function of the physical coordinate.
    /// \param eta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Yd( const double& eta ) const
    {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Yd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime( problem );
    }

    /// Mapping function that provides the second derivative of the
    /// computational Y coordinate as a function of the physical coordinate.
    /// \param eta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Ydd( const double& eta ) const
    {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Ydd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime( problem );
    }
    
    TwoD_Mapped_Node_Mesh()
    {}

    /// ctor
    TwoD_Mapped_Node_Mesh( const double left, const double right, const double bottom, const double top, const std::size_t nx, const std::size_t ny, const std::size_t nvars ) : LEFT( left ), RIGHT(right), BOTTOM(bottom), TOP(top), NX( nx ), NY( ny ), NV( nvars )
    {
      // initialise the storage, but fill these below in init_mapping()
      X = DenseVector<double> (NX,0.0);
      Y = DenseVector<double> (NY,0.0);
      COMP_X = DenseVector<double> (NX,0.0);
      COMP_Y = DenseVector<double> (NY,0.0);
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      VARS = DenseVector<_Type>( NX * NY * NV, 0.0 );
    }

    // Construct the coordinate mapping. The computational mesh is uniform
    // with NX x NY nodes. The physical mesh spans the domain LEFT to RIGHT
    // and BOTTOM to TOP, and its nodes are non-uniformly spaced. The physical
    // nodes are found by inverting the mapping with Newton iteration.
    void init_mapping();

    std::pair<double,double> get_comp_step_sizes() const;

    /// dtor
    virtual ~TwoD_Mapped_Node_Mesh()
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

    /// Access the physical nodal position - as a pair.
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

    /// Get the number of nodes in the two directions of the 2D mesh
    /// \return A pair consisting of the number of nodes in the 2 directions
    std::pair< std::size_t, std::size_t > get_nnodes() const;

    /// Get the number of variables that are stored at each node
    /// \return The number of variables that have data stored at
    /// each nodal point
    std::size_t get_nvars() const;

    /// Access the vector of x-nodal positions
    /// \return A vector of the nodal positions for this mesh
    DenseVector<double>& xnodes();

    /// Access the vector of y-nodal positions
    /// \return A vector of the nodal positions for this mesh
    DenseVector<double>& ynodes();

    /// A simple method for dumping data to a file for gnuplot surface plotting
    /// \param filename The filename to write the data to (will overwrite)
    void dump_gnu( std::string filename ) const;

    /// Normalise all data in the mesh based on one variable.
    /// \param var This var will have its peak (absolute) value as +/-unity
    /// following the normalisation. All other variables will also be rescaled by
    /// the same amount.
    void normalise( const std::size_t& var );

    /// Rescale all values stored in the mapped mesh by a scalar
    /// \param value The scalar that is to multiply all mesh content
    void scale( const _Type& value )
    {
      VARS.scale( value );
    }

    /// Find the maximum stored absolute value in the mesh for a
    /// given variable -- no interpolation is used
    /// \param var The variable index whose maximum is being asked for
    /// \return The value of the maximum (abs value)
    double max( unsigned var )
    {
      double max( 0.0 );
      // step through the nodes
      for ( unsigned nodex = 0; nodex < NX; ++nodex )
      {
        for ( unsigned nodey = 0; nodey < NY; ++nodey )
        {
          if ( std::abs( VARS[ ( nodex * NY + nodey ) * NV + var ] ) > max )
          {
            max = std::abs( VARS[ ( nodex * NY + nodey ) * NV + var ] );
          }
        }
      }
      return max;
    }

  protected:

    // domain boundaries in physical space
    double LEFT,RIGHT,BOTTOM,TOP;
    // we'll store the number of nodes
    std::size_t NX, NY, NV;
    // store computational space (uniform) x nodal points
    // just for convenience since it is uniform we don't need to
    DenseVector<double> COMP_X;
    // store computational space (uniform) y nodal points
    // just for convenience since it is uniform we don't need to
    DenseVector<double> COMP_Y;
    // store physical space (non-uniform) x nodal points
    DenseVector<double> X;
    // store physical space (non-uniform) y nodal points
    DenseVector<double> Y;
    // store the nodal values
    DenseVector<_Type> VARS;
  };

  
  template <typename _Type>
  inline _Type& TwoD_Mapped_Node_Mesh<_Type>::operator()( const std::size_t nodex, const std::size_t nodey, const std::size_t var )
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange( problem, NV, var );
    }
#endif
    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <typename _Type>
  inline const _Type& TwoD_Mapped_Node_Mesh<_Type>::operator()( const std::size_t nodex, const std::size_t nodey, const std::size_t var ) const
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange( problem, NX, nodex, NY, nodey );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange( problem, NV, var );
    }
#endif
    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <typename _Type>
  inline std::pair<double, double> TwoD_Mapped_Node_Mesh<_Type>::coord( const std::size_t nodex, const std::size_t nodey ) const
  {
#ifdef PARANOID
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.coord method is trying to \n";
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

#endif // TWOD_MAPPED_NODE_MESH_H
