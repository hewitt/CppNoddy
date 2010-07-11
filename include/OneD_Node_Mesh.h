/// \file OneD_Node_Mesh.h
/// A specification for a one dimensional mesh object. This is the
/// base class and will be inherited from to form a uniform mesh
/// and non-uniform mesh.

#ifndef ONED_NODE_MESH_H
#define ONED_NODE_MESH_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <DenseVector.h>

namespace CppNoddy
{

  /// A one dimensional mesh utility object. Data can be placed
  /// into the mesh for interpolation to a new mesh or computation
  /// of integrals of the data. The default typing assumes that the
  /// mesh is along the real line. You can store (complex) data across a set
  /// of points in the complex plane using the second typename.
  template < typename _Type, typename _Xtype = double >
  class OneD_Node_Mesh
  {
  public:

    /// Default constructor
    OneD_Node_Mesh()
    {}

    /// ctor for a given nodal distribution
    /// \param nodes The positions of the nodal points
    /// \param nvars The number of variables to store in the mesh
    OneD_Node_Mesh( const DenseVector<_Xtype>& nodes, const std::size_t nvars ) :
        NV( nvars ), X( nodes )
    {
      // set the contents to zero
      VARS = DenseVector<_Type>( NV * X.size(), 0.0 );
    }

    /// Destructor
    virtual ~OneD_Node_Mesh()
    {}

    /// Access the variables at a node as a vector
    /// \param i The index of the node to be accessed
    /// \param var The variable to return the data for
    /// \return The variable 'var' returned for ALL nodal points in the mesh
    _Type& operator()( const std::size_t i, const std::size_t var );

    /// Access the variables at a node as a vector
    /// \param i The index of the node to be accessed
    /// \param var The variable to return the data for
    /// \return The variable 'var' returned for ALL nodal points in the mesh
    const _Type& operator()( const std::size_t i, const std::size_t var ) const;

    /// Const access the nodal position
    /// \param node The nodal position to return
    /// \return The spatial position of this node
    const _Xtype& coord( const std::size_t& node ) const;

    /// Const access the nodal position
    /// \param node The nodal position to return
    /// \return The spatial position of this node
    _Xtype& coord( const std::size_t& node );

    /// Set the variables stored at A SPECIFIED node
    /// \param node The nodal index to be set
    /// \param u The vector of VARIABLES to be written to this nodal point
    void set_nodes_vars( const std::size_t node, const DenseVector<_Type>& u );

    /// Get the variables stored at A SPECIFIED node
    /// \param node The nodal index to be returned
    /// \return The vector of VARIABLES stored at this nodal point
    DenseVector<_Type> get_nodes_vars( const std::size_t& node ) const;

    /// Get the variable data at an interpolated position using
    /// a first order scheme. Doesn't really make sense unless the
    /// data is along the real line.
    /// \param pos The position to interpolate at.
    /// \return A vector of interpolated variables.
    DenseVector<_Type> get_interpolated_vars( const _Type& pos ) const;

    /// \return The number of nodal points in the mesh
    std::size_t get_nnodes() const;

    /// \return The number of variables that have data stored at
    /// each nodal point
    std::size_t get_nvars() const;

    /// \return A vector of the nodal positions for this mesh
    const DenseVector<_Xtype>& nodes() const;

    /// Find a list of approximate locations at which a specified
    /// variable attains a given value.
    /// \param var The variable to be examined for zeros
    /// \param value The value to find
    DenseVector<double> find_roots1( const std::size_t& var, double value = 0.0 ) const;

    /// Integrate over the domain. Typically useful for
    /// finite volume methods.
    /// \param var The variable-index to be integrated over the mesh using a
    /// trapezium rule.
    /// \return The integral value.
    _Type integral2( std::size_t var = 0 ) const;

    /// Compute the integral of the absolute variable squared: |variable|^2.
    /// \param var The variable-index to be integrated
    /// \return The integral of the square of the absolute value.
    _Xtype squared_integral2( std::size_t var = 0 ) const;

    /// Integrate over the domain. Typically useful for
    /// finite volume methods.
    /// \param var The variable-index to be integrated over the mesh using a
    /// trapezium rule.
    /// \return The integral value.
    _Type integral4( std::size_t var = 0 ) const;

    /// For each nodal point we push each variable into a vector
    /// in sequence. So the returned vector has the data
    /// v_00,,...,v_0n,v_10,...,v1n,....v_mn, where the first index
    /// denotes the nodal point 0-m and the second the variable 0-n.
    /// \return A vector of the variables stored in the mesh
    const DenseVector<_Type>& vars_as_vector() const;

    /// Set the variables of this mesh from a vector.
    /// \param vec The vector to be used.
    void set_vars_from_vector( const DenseVector<_Type>& vec );

    /// \return An STL vector of dense vectors of the variables in the mesh
    const std::vector<DenseVector<_Type> >& get_vars() const;

    /// A simple method for dumping data to std::cout
    void dump() const;

    /// Interpolate this mesh data (linearly) into a new
    /// mesh with nodal points defined in the argument list.
    /// \param z The nodal coordinates to be used in the new mesh.
    void remesh1( const DenseVector<_Xtype>& z );

    /// Scale the whole contents of the mesh.
    /// \param x The value to multiply the contents of the mesh by
    void scale( _Type x )
    {
      VARS.scale( x );
    }

    /// A simple method for dumping data to a file for gnuplot
    /// \param filename The filename to write the data to (will overwrite)
    void dump_gnu( std::string filename ) const
    {
      std::ofstream dump;
      dump.open( filename.c_str() );
      dump.precision( 9 );
      dump.setf( std::ios::showpoint );
      dump.setf( std::ios::showpos );
      dump.setf( std::ios::scientific );
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        dump << X[ i ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          dump << VARS[ i * NV + var ] << " ";
        }
        dump << "\n";
      }
    }

  protected:

    // number of variables
    std::size_t NV;
    // store nodal points
    DenseVector<_Xtype> X;
    // store the nodal values
    DenseVector<_Type> VARS;

  };

  // INLINE THE ACCESS METHODS

  template < typename _Type, typename _Xtype >
  inline _Type& OneD_Node_Mesh<_Type, _Xtype>::operator()( const std::size_t i, const std::size_t var )
  {
    return VARS[ i * NV + var ];
  }

  template < typename _Type, typename _Xtype >
  inline const _Type& OneD_Node_Mesh<_Type, _Xtype>::operator()( const std::size_t i, const std::size_t var ) const
  {
    return VARS[ i * NV + var ];
  }

  template < typename _Type, typename _Xtype >
  inline const _Xtype& OneD_Node_Mesh<_Type, _Xtype>::coord( const std::size_t& node ) const
  {
    return X[ node ];
  }

  template < typename _Type, typename _Xtype >
  inline _Xtype& OneD_Node_Mesh<_Type, _Xtype>::coord( const std::size_t& node )
  {
    return X[ node ];
  }
}

#endif // ONED_NODE_MESH_H
