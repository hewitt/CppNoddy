/// \file TwoD_Mapped_Node_Mesh.h
/// A specification for a two dimensional mapped mesh object. Data
/// is stored on a nodal mesh that is non-uniform (in general). This
/// class essentially combines the mapping function for the mesh into
/// the 2D storage of the data.

#ifndef TWOD_MAPPED_NODE_MESH_H
#define TWOD_MAPPED_NODE_MESH_H

#include <vector>
#include <fstream>

#include <OneD_Node_Mesh.h>
#include <DenseVector.h>
#include <Types.h>

namespace CppNoddy {

  /// A two dimensional (mapped) mesh utility object.
  template <typename _Type>
  class TwoD_Mapped_Node_Mesh {
   public:

    /// Mapping function that provides a computational X coordinate
    /// from a physical coordinate.
    /// \param zeta The physical coordinate
    /// \return The corresponding computational coordinate
    virtual double FnComp_X(const double& zeta) const {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_X method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime(problem);
    }

    /// Mapping function that provides the first derivative of the
    /// computational m_X coordinate as a function of the physical coordinate.
    /// \param zeta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Xd(const double& zeta) const {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Xd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime(problem);
    }

    /// Mapping function that provides the second derivative of the
    /// computational X coordinate as a function of the physical coordinate.
    /// \param zeta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Xdd(const double& zeta) const {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Xdd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime(problem);
    }

    /// Mapping function that provides the
    /// computational Y coordinate as a function of the physical coordinate.
    /// \param eta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Y(const double& eta) const {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Y method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime(problem);
    }

    /// Mapping function that provides the first derivative of the
    /// computational Y coordinate as a function of the physical coordinate.
    /// \param eta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Yd(const double& eta) const {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Yd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime(problem);
    }

    /// Mapping function that provides the second derivative of the
    /// computational Y coordinate as a function of the physical coordinate.
    /// \param eta The physical coordinate
    /// \return The corresponding derivative of the computational coordinate
    virtual double FnComp_Ydd(const double& eta) const {
      std::string problem;
      problem = "The TwoD_Mapped_Node_Mesh::FnComp_Ydd method has not been implemented.\n";
      problem += "You have to implement this method to define the mesh.\n";
      throw ExceptionRuntime(problem);
    }

    TwoD_Mapped_Node_Mesh()
    {}

    /// ctor of a blank mesh
    TwoD_Mapped_Node_Mesh(const double left, const double right, const double bottom, const double top, const std::size_t nx, const std::size_t ny, const std::size_t nvars) : m_left(left), m_right(right), m_bottom(bottom), m_top(top), m_nx(nx), m_ny(ny), m_nv(nvars) {
      // initialise the storage, but fill these below in init_mapping()
      m_X = DenseVector<double> (m_nx,0.0);
      m_Y = DenseVector<double> (m_ny,0.0);
      m_compX = DenseVector<double> (m_nx,0.0);
      m_compY = DenseVector<double> (m_ny,0.0);
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      m_vars = DenseVector<_Type>(m_nx * m_ny * m_nv, 0.0);
      // user needs to call init_mapping() to set up m_compX, m_compY after this
    }


    // ctor from a file
    TwoD_Mapped_Node_Mesh(std::string filename, const std::size_t nx, const std::size_t ny, const std::size_t nv)  :
      m_nx(nx), m_ny(ny), m_nv(nv) {
      // need storage for the coordinates
      m_X = DenseVector<double>(m_nx, 0.0);
      m_Y = DenseVector<double>(m_ny, 0.0);
      m_compX = DenseVector<double> (m_nx,0.0);
      m_compY = DenseVector<double> (m_ny,0.0);
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      m_vars = DenseVector<_Type>(m_nx * m_ny * m_nv, 0.0);
      // now read the mesh from the given filename -- this also updates the m_X data
      read(filename, true);
      // set up the private member data on the box size.
      m_left = m_X[0];
      m_right = m_X[m_nx-1];
      m_bottom = m_Y[0];
      m_top = m_Y[m_ny-1];
      // user needs to call init_mapping() to set up m_compX, m_compY after this
    }
    
    
    // Construct the coordinate mapping. The computational mesh is uniform
    // with m_nx x m_ny nodes. The physical mesh spans the domain m_left to m_right
    // and m_bottom to m_top, and its nodes are non-uniformly spaced. The physical
    // nodes are found by inverting the mapping with Newton iteration.
    void init_mapping();

    std::pair<double,double> get_comp_step_sizes() const;

    /// dtor
    virtual ~TwoD_Mapped_Node_Mesh()
    {}

    /// Access operator for a nodal point that returns a vector
    /// \param nodex The nodal index value in the first direction
    /// \param nodey The nodal index value in the second direction
    /// \return The vector of variables stored at the node
    DenseVector<_Type> operator()(const std::size_t nodex, const std::size_t nodey);
    
    /// Access operator for a nodal point/variable in the mesh
    /// \param nodex The nodal index value in the first direction
    /// \param nodey The nodal index value in the second direction
    /// \param var The variable index to be accessed
    _Type& operator()(const std::size_t nodex, const std::size_t nodey, const std::size_t var);

    /// Const access operator for a nodal point/variable in the mesh
    /// \param nodex The nodal index value in the first direction
    /// \param nodey The nodal index value in the second direction
    /// \param var The variable index to be accessed
    const _Type& operator()(const std::size_t nodex, const std::size_t nodey, const std::size_t var) const;

    /// Access the physical nodal position - as a pair.
    /// \param nodex The x nodal position to return
    /// \param nodey The y nodal position to return
    /// \return The spatial position of this node as a pair
    std::pair<double, double> coord(const std::size_t nodex, const std::size_t nodey) const;

    /// Set the variables stored at A SPECIFIED node
    /// \param nodex The x nodal index to be set
    /// \param nodey The y nodal index to be set
    /// \param U The vector of VARIABLES to be written to this nodal point
    void set_nodes_vars(const std::size_t nodex, const std::size_t nodey, const DenseVector<_Type>& U);

    /// Get the variables stored at A SPECIFIED node -- equivalent to mesh(nodex,nodey).
    /// \param nodex The x nodal index to be returned
    /// \param nodey The y nodal index to be returned
    /// \return The vector of VARIABLES stored at this nodal point
    DenseVector<_Type> get_nodes_vars(const std::size_t nodex, const std::size_t nodey) const;

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

    /// A simple method for reading data from a file
    /// \param filename The filename to write the data to (will overwrite)
    /// \param reset Will reset the nodal positions using those from the file
    void read(std::string filename, const bool reset = false);

    /// A simple method for dumping data to a file for gnuplot surface plotting
    /// \param filename The filename to write the data to (will overwrite)
    void dump_gnu(std::string filename) const;

    /// Normalise all data in the mesh based on one variable.
    /// \param var This var will have its peak (absolute) value as +/-unity
    /// following the normalisation. All other variables will also be rescaled by
    /// the same amount.
    void normalise(const std::size_t& var);

    /// Rescale all values stored in the mapped mesh by a scalar
    /// \param value The scalar that is to multiply all mesh content
    void scale(const _Type& value) {
      m_vars.scale(value);
    }

    /// Find the maximum stored absolute value in the mesh for a
    /// given variable -- no interpolation is used
    /// \param var The variable index whose maximum is being asked for
    /// \return The value of the maximum (abs value)
    double max(unsigned var) {
      double max(0.0);
      // step through the nodes
      for(unsigned nodex = 0; nodex < m_nx; ++nodex) {
        for(unsigned nodey = 0; nodey < m_ny; ++nodey) {
          if(std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ]) > max) {
            max = std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ]);
          }
        }
      }
      return max;
    }
    
    /// Get a bilinearly interpolated value at a specified point
    /// \param x Physical (unmapped) x-coordinate in the 2D mesh
    /// \param y Physical (unmapped) y-coordinate in the 2D mesh
    /// \return A vector of bilinearly interpolated values
    DenseVector<_Type> get_interpolated_vars(const double& x, const double& y) 
    {
      double compX = FnComp_X(x);
      double compY = FnComp_Y(y);
      const double tol(1.e-8);
      // check start and end
      if ((compX < m_compX[0] - tol) || (compX > m_compX[m_nx-1] + tol)) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method has been called with \n";
        problem += " an x coordinate that lies outside the mesh. \n";
        throw ExceptionRuntime(problem);
      }
      // check start and end
      if ((compY < m_compY[0] - tol) || (compY > m_compY[m_ny-1] + tol)) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method has been called with \n";
        problem += " a y coordinate that lies outside the mesh. \n";
        throw ExceptionRuntime(problem);
      }
      int bottom_j(-1);
      for (unsigned j = 0; j < m_ny-1; ++j) {
        if ((compY >= m_compY[j] - tol) && (compY <= m_compY[j+1] + tol)) {
          bottom_j = j;
        }
      }
      if (bottom_j == -1) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method is broken.\n";
        throw ExceptionRuntime(problem);
      }
      //
      OneD_Node_Mesh<_Type> bottom_row( m_compX, m_nv );
      OneD_Node_Mesh<_Type> top_row( m_compX, m_nv );
      for (std::size_t i = 0; i < m_nx; ++i ) {
        bottom_row.set_nodes_vars(i, this -> get_nodes_vars(i, bottom_j));
        top_row.set_nodes_vars(i, this -> get_nodes_vars(i, bottom_j+1));
      }
      const double compY1 = m_compY[ bottom_j ];
      const double compY2 = m_compY[ bottom_j+1 ];
      DenseVector<_Type> result =
        top_row.get_interpolated_vars(compX)*(compY-compY1)/(compY2-compY1)
        + bottom_row.get_interpolated_vars(compX)*(compY2-compY)/(compY2-compY1);
      return result;
    }

    

   protected:

    // domain boundaries in physical space
    double m_left,m_right,m_bottom,m_top;
    // we'll store the number of nodes
    std::size_t m_nx, m_ny, m_nv;
    // store computational space (uniform) x nodal points
    // just for convenience since it is uniform we don't need to
    DenseVector<double> m_compX;
    // store computational space (uniform) y nodal points
    // just for convenience since it is uniform we don't need to
    DenseVector<double> m_compY;
    // store physical space (non-uniform) x nodal points
    DenseVector<double> m_X;
    // store physical space (non-uniform) y nodal points
    DenseVector<double> m_Y;
    // store the nodal values
    DenseVector<_Type> m_vars;
  };

  template <typename _Type>
  inline DenseVector<_Type> TwoD_Mapped_Node_Mesh<_Type>::operator()(const std::size_t nodex, const std::size_t nodey ) {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
#endif
    return get_nodes_vars(nodex,nodey);
  }

  
  template <typename _Type>
  inline _Type& TwoD_Mapped_Node_Mesh<_Type>::operator()(const std::size_t nodex, const std::size_t nodey, const std::size_t var) {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
    if(var > m_nv - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nv, var);
    }
#endif
    return m_vars[(nodex * m_ny + nodey) * m_nv + var ];
  }

  template <typename _Type>
  inline const _Type& TwoD_Mapped_Node_Mesh<_Type>::operator()(const std::size_t nodex, const std::size_t nodey, const std::size_t var) const {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
    if(var > m_nv - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nv, var);
    }
#endif
    return m_vars[(nodex * m_ny + nodey) * m_nv + var ];
  }

  template <typename _Type>
  inline std::pair<double, double> TwoD_Mapped_Node_Mesh<_Type>::coord(const std::size_t nodex, const std::size_t nodey) const {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.coord method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
#endif
    std::pair< double, double > pos;
    pos.first = m_X[ nodex ];
    pos.second = m_Y[ nodey ];
    return pos;
  }

}

#endif // TWOD_MAPPED_NODE_MESH_H
