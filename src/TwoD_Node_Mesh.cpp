/// \file TwoD_Node_Mesh.cpp
/// Implementation of a two dimensional mesh object. Data
/// is stored on a nodal mesh.

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Exceptions.h>
#include <DenseVector.h>
#include <DenseMatrix.h>
#include <TwoD_Node_Mesh.h>
#include <Utility.h>

namespace CppNoddy {

  template <typename _Type>
  void TwoD_Node_Mesh<_Type>::set_nodes_vars(const std::size_t nodex, const std::size_t nodey, const DenseVector<_Type>& U) {
#ifdef PARANOID
    if(U.size() > m_nv) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.set_nodes_vars method is trying to use a \n";
      problem += " vector that has more entries than variables stored in the mesh. \n";
      throw ExceptionRuntime(problem);
    }
#endif
    // assign contents of U to the member data
    std::size_t offset((nodex * m_ny + nodey) * m_nv);
    for(std::size_t var = 0; var < m_nv; ++var) {
      m_vars[ offset++ ] = U[ var ];
    }
  }

  template <typename _Type>
  DenseVector<_Type> TwoD_Node_Mesh<_Type>::get_nodes_vars(const std::size_t nodex, const std::size_t nodey) const {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.get_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
#endif
    // construct a vector with m_nv elements starting from a pointer
    DenseVector<_Type> nodes_vars(m_nv, &m_vars[(nodex * m_ny + nodey) * m_nv ]);
    return nodes_vars;
  }

  template <typename _Type>
  void TwoD_Node_Mesh<_Type>::assign(const _Type elt) {
    m_vars.assign(m_nx * m_ny * m_nv, elt);
  }

  template <typename _Type>
  std::pair< std::size_t, std::size_t > TwoD_Node_Mesh<_Type>::get_nnodes() const {
    std::pair< std::size_t, std::size_t > nodes;
    nodes.first = m_nx;
    nodes.second = m_ny;
    return nodes;
  }

  template <typename _Type>
  std::size_t TwoD_Node_Mesh<_Type>::get_nvars() const {
    return m_nv;
  }

  template <typename _Type>
  const DenseVector<double>& TwoD_Node_Mesh<_Type>::xnodes() const {
    return m_X;
  }

  template <typename _Type>
  const DenseVector<double>& TwoD_Node_Mesh<_Type>::ynodes() const {
    return m_Y;
  }

  template <typename _Type>
  DenseMatrix<_Type> TwoD_Node_Mesh<_Type>::get_var_as_matrix(std::size_t var) const {
#ifdef PARANOID
    if(var > m_nv - 1) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.get_var_as_matrix method is trying to use a \n";
      problem += " variable index bigger than the number of variables in the mesh. \n";
      throw ExceptionRange(problem, m_nv, var);
    }
#endif
    DenseMatrix<_Type> temp(m_nx, m_ny, 0.0);
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        temp(i, j) = m_vars[(i * m_ny + j) * m_nv + var ];
      }
    }
    return temp;
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::remesh1(const DenseVector<double>& newX, const DenseVector<double>& newY) {
#ifdef PARANOID
    // check start & end
    if(std::abs(m_X[ 0 ] - newX[ 0 ]) > 1.e-10 ||
        std::abs(m_X[ m_X.size() - 1 ] - newX[ newX.size() - 1 ]) > 1.e-10) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.remesh1 method has been called with \n";
      problem += " a passed X coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw ExceptionRuntime(problem);
    }
    // check monotonic node positions
    for(std::size_t i = 0; i < newX.size() - 1; ++i) {
      if(newX[ i ] >= newX[ i + 1 ]) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.remesh1 method has been passed \n";
        problem += " a non-monotonic X coordinate vector. \n";
        problem += Utility::stringify(newX[ i ], 6) + " vs. " + Utility::stringify(newX[ i + 1 ], 6);
        throw ExceptionRuntime(problem);
      }
    }
    // check start and end
    if(std::abs(m_Y[ 0 ] - newY[ 0 ]) > 1.e-10 ||
        std::abs(m_Y[ m_Y.size() - 1 ] - newY[ newY.size() - 1 ]) > 1.e-10) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.remesh1 method has been called with \n";
      problem += " a passed Y coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw ExceptionRuntime(problem);
    }
    // check monotonic node positions
    for(std::size_t i = 0; i < newY.size() - 1; ++i) {
      if(newY[ i ] >= newY[ i + 1 ]) {
        std::string problem;
        problem = " The TwoD_Node_Mesh.remesh1 method has been passed \n";
        problem += " a non-monotonic Y coordinate vector. \n";
        problem += Utility::stringify(newY[ i ], 6) + " vs. " + Utility::stringify(newY[ i + 1 ], 6);
        throw ExceptionRuntime(problem);
      }
    }
#endif

    // new variables storage
    DenseVector<_Type> newvars(newX.size() * newY.size() * m_nv, 0.0);

    // left boundary
    {
      std::size_t xnode(0);
      // bottom left corner copy
      for(unsigned var = 0; var < m_nv; ++var) {
        newvars[(xnode * newY.size() + 0) * m_nv + var ] = get_nodes_vars(0, 0)[ var ];
      }
      for(std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode) {
        std::size_t left_i(0);    // bracketing index
        std::size_t below_j(0);   // bracketing index
        double deltaY(0.0);
        // loop through the source mesh and find the bracket-nodes
        for(std::size_t j = 0; j < m_Y.size() - 1; ++j) {
          if((m_Y[ j ] <= newY[ ynode ]) && (newY[ ynode ] < m_Y[ j + 1 ])) {
            below_j = j;
            deltaY = newY[ ynode ] - m_Y[ j ];
          }
        }
        DenseVector<_Type> dvarsdY = (get_nodes_vars(left_i, below_j + 1) - get_nodes_vars(left_i, below_j))
                                     / (coord(left_i, below_j + 1).second - coord(left_i, below_j).second);
        DenseVector<_Type> interpolated_vars = get_nodes_vars(left_i, below_j) + dvarsdY * deltaY;
        for(unsigned var = 0; var < m_nv; ++var) {
          newvars[(xnode * newY.size() + ynode) * m_nv + var ] = interpolated_vars[ var ];
        }
      }
      // top left corner copy
      for(unsigned var = 0; var < m_nv; ++var) {
        newvars[(xnode * newY.size() + newY.size() - 1) * m_nv + var ] = get_nodes_vars(0, m_ny - 1)[ var ];
      }
    }
    // right boundary
    {
      std::size_t xnode(newX.size() - 1);
      // bottom right corner copy
      for(unsigned var = 0; var < m_nv; ++var) {
        newvars[(xnode * newY.size() + 0) * m_nv + var ] = get_nodes_vars(m_nx - 1, 0)[ var ];
      }
      for(std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode) {
        std::size_t left_i(m_X.size() - 1);    // bracketing index
        std::size_t below_j(0);   // bracketing index
        double deltaY(0.0);
        // loop through the source mesh and find the bracket-nodes
        for(std::size_t j = 0; j < m_Y.size() - 1; ++j) {
          if((m_Y[ j ] <= newY[ ynode ]) && (newY[ ynode ] < m_Y[ j + 1 ])) {
            below_j = j;
            deltaY = newY[ ynode ] - m_Y[ j ];
          }
        }
        DenseVector<_Type> dvarsdY = (get_nodes_vars(left_i, below_j + 1) - get_nodes_vars(left_i, below_j))
                                     / (coord(left_i, below_j + 1).second - coord(left_i, below_j).second);
        DenseVector<_Type> interpolated_vars = get_nodes_vars(left_i, below_j) + dvarsdY * deltaY;
        for(unsigned var = 0; var < m_nv; ++var) {
          newvars[(xnode * newY.size() + ynode) * m_nv + var ] = interpolated_vars[ var ];
        }
      }
      // bottom right corner copy
      for(unsigned var = 0; var < m_nv; ++var) {
        newvars[(xnode * newY.size() + newY.size() - 1) * m_nv + var ] = get_nodes_vars(m_nx - 1, m_ny - 1)[ var ];
      }
    }
    // bottom boundary
    {
      std::size_t ynode(0);
      for(std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode) {
        std::size_t left_i(0);    // bracketing index
        std::size_t below_j(0);   // bracketing index
        double deltaX(0.0);
        // loop through the source mesh and find the bracket-nodes
        for(std::size_t i = 0; i < m_X.size() - 1; ++i) {
          if((m_X[ i ] <= newX[ xnode ]) && (newX[ xnode ] < m_X[ i + 1 ])) {
            left_i = i;
            deltaX = newX[ xnode ] - m_X[ i ];
          }
        }
        DenseVector<_Type> dvarsdX = (get_nodes_vars(left_i + 1, below_j) - get_nodes_vars(left_i, below_j))
                                     / (coord(left_i + 1, below_j).first - coord(left_i, below_j).first);
        DenseVector<_Type> interpolated_vars = get_nodes_vars(left_i, below_j) + dvarsdX * deltaX;
        for(unsigned var = 0; var < m_nv; ++var) {
          newvars[(xnode * newY.size() + ynode) * m_nv + var ] = interpolated_vars[ var ];
        }
      }
    }
    // top boundary
    {
      std::size_t ynode(newY.size() - 1);
      for(std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode) {
        std::size_t left_i(0);    // bracketing index
        std::size_t below_j(m_Y.size() - 1);   // bracketing index
        double deltaX(0.0);
        // loop through the source mesh and find the bracket-nodes
        for(std::size_t i = 0; i < m_X.size() - 1; ++i) {
          if((m_X[ i ] <= newX[ xnode ]) && (newX[ xnode ] < m_X[ i + 1 ])) {
            left_i = i;
            deltaX = newX[ xnode ] - m_X[ i ];
          }
        }
        DenseVector<_Type> dvarsdX = (get_nodes_vars(left_i + 1, below_j) - get_nodes_vars(left_i, below_j))
                                     / (coord(left_i + 1, below_j).first - coord(left_i, below_j).first);
        DenseVector<_Type> interpolated_vars = get_nodes_vars(left_i, below_j) + dvarsdX * deltaX;
        for(unsigned var = 0; var < m_nv; ++var) {
          newvars[(xnode * newY.size() + ynode) * m_nv + var ] = interpolated_vars[ var ];
        }
      }
    }
    // loop thru interior nodes of the destination mesh one node at a time
    for(std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode) {
      for(std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode) {
        std::size_t left_i(0);    // bracketing index
        std::size_t below_j(0);   // bracketing index
        // loop through the source mesh and find the bracket-nodes
        for(std::size_t i = 0; i < m_X.size() - 1; ++i) {
          if((m_X[ i ] <= newX[ xnode ]) && (newX[ xnode ] < m_X[ i + 1 ])) {
            left_i = i;
          }
        }
        // loop through the source mesh and find the bracket-nodes
        for(std::size_t j = 0; j < m_Y.size() - 1; ++j) {
          if((m_Y[ j ] <= newY[ ynode ]) && (newY[ ynode ] < m_Y[ j + 1 ])) {
            below_j = j;
          }
        }
        DenseVector<_Type> dvarsdX = (get_nodes_vars(left_i + 1, below_j) - get_nodes_vars(left_i, below_j))
                                     / (coord(left_i + 1, below_j).first - coord(left_i, below_j).first);
        DenseVector<_Type> dvarsdY = (get_nodes_vars(left_i, below_j + 1) - get_nodes_vars(left_i, below_j))
                                     / (coord(left_i, below_j + 1).second - coord(left_i, below_j).second);

        DenseVector<_Type> interpolated_vars_bottom =
          (get_nodes_vars(left_i, below_j) * (coord(left_i + 1, below_j).first - newX[ xnode ])
           + get_nodes_vars(left_i + 1, below_j) * (newX[ xnode ] - coord(left_i, below_j).first)) /
          (coord(left_i + 1, below_j).first - coord(left_i, below_j).first);

        DenseVector<_Type> interpolated_vars_top =
          (get_nodes_vars(left_i, below_j + 1) * (coord(left_i + 1, below_j + 1).first - newX[ xnode ])
           + get_nodes_vars(left_i + 1, below_j + 1) * (newX[ xnode ] - coord(left_i, below_j + 1).first)) /
          (coord(left_i + 1, below_j + 1).first - coord(left_i, below_j + 1).first);

        DenseVector<_Type> interpolated_vars =
          (interpolated_vars_bottom * (coord(left_i, below_j + 1).second - newY[ ynode ])
           +  interpolated_vars_top * (newY[ ynode ] - coord(left_i, below_j).second)) /
          (coord(left_i, below_j + 1).second - coord(left_i, below_j).second);

        for(unsigned var = 0; var < m_nv; ++var) {
          newvars[(xnode * newY.size() + ynode) * m_nv + var ] = interpolated_vars[ var ];
        }
      }
    }
    // finally replace the old nodes with the new ones
    m_X = newX;
    m_Y = newY;
    m_nx = newX.size();
    m_ny = newY.size();
    m_vars = newvars;
  }

  template<typename _Type>
  OneD_Node_Mesh<_Type> TwoD_Node_Mesh<_Type>::get_xsection_at_xnode(const std::size_t nodex) const {
    OneD_Node_Mesh<_Type> xsection(m_Y, m_nv);
    for(std::size_t nodey = 0; nodey < m_ny; ++nodey) {
      xsection.set_nodes_vars(nodey, this -> get_nodes_vars(nodex, nodey));
    }
    return xsection;
  }

  template<typename _Type>
  OneD_Node_Mesh<_Type> TwoD_Node_Mesh<_Type>::get_xsection_at_ynode(const std::size_t nodey) const {
    OneD_Node_Mesh<_Type> xsection(m_X, m_nv);
    for(std::size_t nodex = 0; nodex < m_nx; ++nodex) {
      xsection.set_nodes_vars(nodex, this -> get_nodes_vars(nodex, nodey));
    }
    return xsection;
  }


  template <typename _Type>
  void TwoD_Node_Mesh<_Type>::dump() const {
    for(std::size_t var = 0; var < m_nv; ++var) {
      std::cout << "Variable : " << var << "\n";
      std::cout << " x = ";
      for(std::size_t i = 0; i < m_nx; ++i) {
        std::cout << m_X[ i ] << ", ";
      }

      std::cout << "\n";
      for(std::size_t j = 0; j < m_ny; ++j) {
        std::cout << " y = " << m_Y[ j ] << "\n";
        for(std::size_t i = 0; i < m_nx; ++i) {
          std::cout << m_vars[(i * m_ny + j) * m_nv + var ] << ", ";
        }
        std::cout << "\n";
      }
    }
  }


  template<>
  void TwoD_Node_Mesh<double>::normalise(const std::size_t& var) {
    double maxval(max_abs(var));
    m_vars.scale(1./maxval);
  }

  template<>
  void TwoD_Node_Mesh<D_complex>::normalise(const std::size_t& var) {
    unsigned max_nx(0);
    unsigned max_ny(0);
    double max(0.0);
    // step through the nodes
    for(unsigned nodex = 0; nodex < m_X.size(); ++nodex) {
      for(unsigned nodey = 0; nodey < m_Y.size(); ++nodey) {
        if(std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ]) > max) {
          max = std::abs(m_vars[(nodex * m_ny + nodey) * m_nv + var ]);
          max_nx = nodex;
          max_ny = nodey;
        }
      }
    }
    D_complex factor(m_vars[(max_nx * m_ny + max_ny) * m_nv + var ]);
    std::cout << "[DEBUG] Normalise factor = " << factor << "\n";
    m_vars.scale(1./factor);
    // D_complex max_elt = m_vars[(max_nx * m_ny + max_ny) * m_nv + var ];
    // std::cout << "[DEBUG] Normalise max_elt = " << max_elt << "\n";
    // // to here will give |variable|=1 but in general v_r and v+i .ne. 0
    // // additionally we can scale by a factor of unit magnitude
    // // in order to make v_r=1 and v_i=0.
    // m_vars.scale( std::conj(max_elt) );
    // std::cout << "[DEBUG] Normalise max_elt = " << m_vars[(max_nx * m_ny + max_ny) * m_nv + var ] << "\n";
  }

  // specialised because obviously a double mesh element won't be able
  // to call the .real() method.
  template<>
  double TwoD_Node_Mesh<D_complex>::max_abs_real_part(unsigned var) {
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

  template<>
  double TwoD_Node_Mesh<double>::max_abs_real_part(unsigned var) {
#ifdef PARANOID
    // check start & end
    std::string problem;
    problem = " The TwoD_Node_Mesh.max_abs_real_part method has been called for \n";
    problem += " a mesh with element type of double (rather than D_complex).";
    throw ExceptionRuntime(problem);
#endif
    // just call the max_abs
    return max_abs(var);
  }

  template<>
  void TwoD_Node_Mesh<double>::dump_gnu(std::string filename) const {
    std::ofstream dump;
    dump.open(filename.c_str());
    dump.precision(15);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    //
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        dump << m_X[ i ] << " " << m_Y[ j ] << " ";
        for(std::size_t var = 0; var < m_nv; ++var) {
          dump << m_vars[(i * m_ny + j) * m_nv + var ] << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

  template <>
  void TwoD_Node_Mesh<D_complex>::dump_gnu(std::string filename) const {
    std::ofstream dump;
    dump.open(filename.c_str());
    dump.precision(15);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    //
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        dump << m_X[ i ] << " " << m_Y[ j ] << " ";
        for(std::size_t var = 0; var < m_nv; ++var) {
          dump << real(m_vars[(i * m_ny + j) * m_nv + var ]) << " ";
          dump << imag(m_vars[(i * m_ny + j) * m_nv + var ]) << " ";
        }
        dump << "\n";
      }
      dump << "\n";
    }
    dump.close();
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::dump_var(std::string filename, const unsigned var) const {
    std::ofstream dump;
    dump.open(filename.c_str());
    dump.precision(15);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        dump << m_vars[(i * m_ny + j) * m_nv + var ] << "\n";
      }
    }
  }

  template< typename _Type>
  void TwoD_Node_Mesh<_Type>::dump(std::string filename) const {
    std::ofstream dump;
    dump.open(filename.c_str());
    dump.precision(15);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    //dump << m_nx << " " << m_ny << " " << m_nv << "\n";
    dump.precision(9);
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        dump << m_X[ i ] << " " << m_Y[ j ] << " ";
        for(std::size_t var = 0; var < m_nv; ++var) {
          dump << m_vars[(i * m_ny + j) * m_nv + var ] << " ";
        }
        dump << "\n";
      }
    }
  }

  template<>
  void TwoD_Node_Mesh<double>::read(std::string filename, bool reset) {
    std::ifstream dump;
    dump.open(filename.c_str());
    if(dump.good() != true) {
      std::string problem;
      problem = " The TwoD_Node_Mesh.read method is trying to read a \n";
      problem += " file (" + filename + ") that doesn't exist.\n";
      throw ExceptionRuntime(problem);
    }
    dump.precision(15);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        double x, y;
        dump >> x;
        dump >> y;
        for(std::size_t var = 0; var < m_nv; ++var) {
          double value;
          dump >> value;
          m_vars[(i * m_ny + j) * m_nv + var ] = value;
        }
        if(reset != true) {
          // if not reseting the mesh we should check the node positions
          if((std::fabs(x - m_X[ i ]) > 1.e-6) || (std::fabs(y - m_Y[ j ]) > 1.e-6)) {
            std::cout << " Read x = " << x << " Expected x = " << m_X[ i ] << "; Read y = " << y << " Expected y = " << m_Y[ j ] << " \n";
            std::cout << " Absolute differences are " << fabs(x - m_X[i]) << " and " << fabs(y - m_Y[j]) << "\n";
            std::string problem;
            problem = " The TwoD_Node_Mesh.read method is trying to read a \n";
            problem += " file whose nodal points are in a different position. \n";
            throw ExceptionRuntime(problem);
          }
        } else {
          m_X[ i ] = x;
          m_Y[ j ] = y;
        }
      }
    }
  }


  template<>
  void TwoD_Node_Mesh<D_complex>::read(std::string filename, bool reset) {
    std::ifstream dump;
    dump.open(filename.c_str());
    dump.precision(15);
    dump.setf(std::ios::showpoint);
    dump.setf(std::ios::showpos);
    dump.setf(std::ios::scientific);
    //
    // 18/06/2017: switched i and j below for consistency with double
    //
    for(std::size_t i = 0; i < m_nx; ++i) {
      for(std::size_t j = 0; j < m_ny; ++j) {
        double x, y;
        dump >> x;
        dump >> y;
        for(std::size_t var = 0; var < m_nv; ++var) {
          double value_r, value_i;
          dump >> value_r;
          dump >> value_i;
          m_vars[(i * m_ny + j) * m_nv + var ] = D_complex(value_r, value_i);
        }
        if(reset != true) {
          // if not reseting the mesh we should check the node positions
          if((std::fabs(x - m_X[ i ]) > 1.e-6) || (std::fabs(y - m_Y[ j ]) > 1.e-6)) {
            std::cout << " Read x = " << x << " Expected x = " << m_X[ i ] << "; Read y = " << y << " Expected y = " << m_Y[ j ] << " \n";
            std::cout << " Absolute differences are " << fabs(x - m_X[i]) << " and " << fabs(y - m_Y[j]) << "\n";
            std::string problem;
            problem = " The TwoD_Node_Mesh.read method is trying to read a \n";
            problem += " file whose nodal points are in a different position. \n";
            throw ExceptionRuntime(problem);
          }
        } else {
          m_X[ i ] = x;
          m_Y[ j ] = y;
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
