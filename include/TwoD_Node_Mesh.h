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
#include <Types.h>

namespace CppNoddy {

  /// A two dimensional mesh utility object.
  template <typename _Type>
  class TwoD_Node_Mesh {
   public:

    TwoD_Node_Mesh()
    {}

    /// ctor
    TwoD_Node_Mesh(const DenseVector<double>& x_nodes, const DenseVector<double>& y_nodes,
                   const std::size_t nvars) :
      m_nx(x_nodes.size()), m_ny(y_nodes.size()), m_nv(nvars), m_X(x_nodes), m_Y(y_nodes) {
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      m_vars = DenseVector<_Type>(m_nx * m_ny * m_nv, 0.0);
    }

    /// ctor
    TwoD_Node_Mesh(const double left, const double right, const double bottom, const double top,
                   const std::size_t nx, const std::size_t ny, const std::size_t nvars) :
      m_nx(nx), m_ny(ny), m_nv(nvars) {
      {
        m_X.reserve(m_nx);
        const double delta = (right - left) / (m_nx - 1);
        for(std::size_t i = 0; i < m_nx; ++i) {
          m_X.push_back(left + delta * i);
        }
      }
      {
        m_Y.reserve(m_ny);
        const double delta = (top - bottom) / (m_ny - 1);
        for(std::size_t i = 0; i < m_ny; ++i) {
          m_Y.push_back(bottom + delta * i);
        }
      }
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      m_vars = DenseVector<_Type>(m_nx * m_ny * m_nv, 0.0);
    }

    // ctor from a file
    TwoD_Node_Mesh(std::string filename, const std::size_t nx, const std::size_t ny, const std::size_t nv)  :
      m_nx(nx), m_ny(ny), m_nv(nv) {
      // need storage for the coordinates
      m_X = DenseVector<double>(nx, 0.0);
      m_Y = DenseVector<double>(ny, 0.0);
      // we'll store the data as ( x, y, v ) ->  x * ny * nv + y * nv + v
      m_vars = DenseVector<_Type>(nx * ny * nv, 0.0);
      // now read the mesh from the given filename
      read(filename, true);
    }

    /// dtor
    virtual ~TwoD_Node_Mesh()
    {}

    /// Access operator for a nodal point that returns a vector
    /// \param nodex The nodal index value in the first direction
    /// \param nodey The nodal index value in the second direction
    /// \return The vector of variables stored at the node
    DenseVector<_Type> operator()(const std::size_t nodex, const std::size_t nodey );
    
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

    /// Access the nodal position - as a pair.
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

    /// Get a cross section of the 2D mesh at a specified (constant) x node
    /// \param nodex The x nodal index at which the cross section is to be taken
    /// \return A 1D nodal mesh
    OneD_Node_Mesh<_Type> get_xsection_at_xnode(const std::size_t nodex) const;

    /// Get a cross section of the 2D mesh at a specified (constant) y node
    /// \param nodey The y nodal index at which the cross section is to be taken
    /// \return A 1D nodal mesh
    OneD_Node_Mesh<_Type> get_xsection_at_ynode(const std::size_t nodey) const;

    /// Get a cross section of the 2D mesh at a specified (constant) x node
    /// \param nodex The x nodal index at which the cross section is to be taken
    /// \return A 1D nodal mesh
    OneD_Node_Mesh<_Type> get_xsection_at_x1(const double x) const {
      unsigned I(0);
      OneD_Node_Mesh<_Type> xsection(m_Y, m_nv);
      for(unsigned i = 0; i<m_nx-1; ++i) {
        if((m_X[i]< x) && (m_X[i+1]>x)) {
          I=i;
        }
      }
      double dx_ratio((x-m_X[I])/(m_X[I+1]-m_X[I]));
      for(unsigned j = 0; j<m_ny; ++j) {
        for(unsigned var = 0; var<m_nv; ++var) {
          xsection(j,var) = this->operator()(I,j,var)+(this->operator()(I+1,j,var)-this->operator()(I,j,var))*dx_ratio;
        }
      }
      return xsection;
    }

    /// Assign an element to all entries in the mesh
    /// \param elt The element to be assigned to the mesh
    void assign(const _Type elt);

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

    /// Return a matrix corresponding to each nodal point in the mesh
    /// Each matrix element will contain a specified variable number
    /// \param var The variable number to be accessed
    /// \return A dense matrix of the specified variable
    DenseMatrix<_Type> get_var_as_matrix(std::size_t var) const;

    /// Interpolate this mesh data (bilinearly) into a new
    /// mesh with nodal points defined in the argument list.
    /// Not written to be efficient, so you probably don't want
    /// to do any repeated calls with this method.
    /// \param newX The x-nodal coordinates to be used in the new mesh.
    /// \param newY The y-nodal coordinates to be used in the new mesh.
    void remesh1(const DenseVector<double>& newX, const DenseVector<double>& newY);

    /// A simple method for dumping data to std::cout
    void dump() const;

    /// A simple method for dumping data to a file
    /// \param filename The filename to write the data to (will overwrite)
    void dump(std::string filename) const;

    /// A simple method for dumping a single variable to a file with no nodal information
    /// \param filename The filename to write the data to (will overwrite)
    /// \param var The index of the variable to be dumped to output
    void dump_var(std::string filename, const unsigned var) const;

    /// A simple method for reading data from a file
    /// \param filename The filename to write the data to (will overwrite)
    /// \param reset Will reset the nodal positions using those from the file
    void read(std::string filename, const bool reset = false);

    /// A simple method for dumping data to a file for gnuplot surface plotting
    /// \param filename The filename to write the data to (will overwrite)
    void dump_gnu(std::string filename) const;

    /// Normalise all data in the mesh based on one variable.
    /// \param var This var will have its peak (absolute) value as +/-unity following
    /// the normalisation. All other variables will also be rescaled by
    /// the same amount.
    void normalise(const std::size_t& var);

    void scale(const _Type& value) {
      m_vars.scale(value);
    }

    void normalise_real_part(const std::size_t& var) {
      double maxval(max_real_part(var));
      m_vars.scale(1./maxval);
    }

    /// Find the maximum stored absolute value in the mesh for a given variable -- no interpolation is used
    /// \param var The variable index whose maximum is being asked for
    /// \return The value of the maximum (abs value)
    double max_real_part(unsigned var) {
      double max(0.0);
      // step through the nodes
      for(unsigned nodex = 0; nodex < m_X.size(); ++nodex) {
        for(unsigned nodey = 0; nodey < m_Y.size(); ++nodey) {
          if(std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ].real()) > max) {
            max = std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ].real());
          }
        }
      }
      return max;
    }

    /// Find the maximum stored absolute value in the mesh for a given variable -- no interpolation is used
    /// \param var The variable index whose maximum is being asked for
    /// \return The value of the maximum (abs value)
    double max(unsigned var) {
      double max(0.0);
      // step through the nodes
      for(unsigned nodex = 0; nodex < m_X.size(); ++nodex) {
        for(unsigned nodey = 0; nodey < m_Y.size(); ++nodey) {
          if(std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ]) > max) {
            max = std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ]);
          }
        }
      }
      return max;
    }


    /// Get a bilinearly interpolated value at a specified point
    /// \param x x-coordinate in the 2D mesh
    /// \param y y-coordinate in the 2D mesh
    /// \return A vector of bilinearly interpolated values
    DenseVector<_Type> get_interpolated_vars(const double& x, const double& y) {
      const double tol(1.e-10);
      // check start and end
      if((x < m_X[0] - tol) || (x > m_X[m_nx-1] + tol)) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method has been called with \n";
        problem += " an x coordinate that lies outside the mesh. \n";
        throw ExceptionRuntime(problem);
      }
      // check start and end
      if((y < m_Y[0] - tol) || (y > m_Y[m_ny-1] + tol)) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method has been called with \n";
        problem += " a y coordinate that lies outside the mesh. \n";
        throw ExceptionRuntime(problem);
      }
      int bottom_j(-1);
      for(unsigned j = 0; j < m_ny-1; ++j) {
        if((y >= m_Y[j] - tol) && (y <= m_Y[j+1] + tol)) {
          bottom_j = j;
        }
        //if ( abs(y-Y[j]) < tol )
        //{
        //bottom_j = j;
        //}
        //if ( abs(y-m_Y[j+1]) < tol )
        //{
        //bottom_j = j+1;
        //}
      }
      //std::cout << y << " " << m_Y[bottom_j] << " " << m_Y[bottom_j+1] << "\n";
      if(bottom_j == -1) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.get_interpolated_vars method is broken.\n";
        throw ExceptionRuntime(problem);
      }
      //
      OneD_Node_Mesh<_Type> bottom_row = get_xsection_at_ynode(bottom_j);
      OneD_Node_Mesh<_Type> top_row = get_xsection_at_ynode(bottom_j+1);
      const double y1 = m_Y[ bottom_j ];
      const double y2 = m_Y[ bottom_j+1 ];
      DenseVector<_Type> result = top_row.get_interpolated_vars(x)*(y-y1)/(y2-y1)
                                  + bottom_row.get_interpolated_vars(x)*(y2-y)/(y2-y1);
      //std::cout << "x,y,interp: " << x << " " << y << " " << result[0] << "\n";
      return result;
    }


   protected:

    // we'll store the number of nodes
    std::size_t m_nx, m_ny, m_nv;
    // store x nodal points
    DenseVector<double> m_X;
    // store y nodal points
    DenseVector<double> m_Y;
    // store the nodal values
    DenseVector<_Type> m_vars;
  };


  template <typename _Type>
  inline DenseVector<_Type> TwoD_Node_Mesh<_Type>::operator()(const std::size_t nodex, const std::size_t nodey ) {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
#endif
    return get_nodes_vars(nodex, nodey);
  }


  template <typename _Type>
  inline _Type& TwoD_Node_Mesh<_Type>::operator()(const std::size_t nodex, const std::size_t nodey, const std::size_t var) {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
    if(var > m_nv - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nv, var);
    }
#endif
    return m_vars[(nodex * m_ny + nodey) * m_nv + var ];
  }

  template <typename _Type>
  inline const _Type& TwoD_Node_Mesh<_Type>::operator()(const std::size_t nodex, const std::size_t nodey, const std::size_t var) const {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
    if(var > m_nv - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nv, var);
    }
#endif
    return m_vars[(nodex * m_ny + nodey) * m_nv + var ];
  }

  template <typename _Type>
  inline std::pair<double, double> TwoD_Node_Mesh<_Type>::coord(const std::size_t nodex, const std::size_t nodey) const {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.coord method is trying to \n";
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

#endif // TWOD_NODE_MESH_H
