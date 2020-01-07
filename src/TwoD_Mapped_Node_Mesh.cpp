/// \file TwoD_Mapped_Node_Mesh.cpp
/// Implementation of a two dimensional (mapped) mesh object. Data
/// is stored on a (mapped) nodal mesh.

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <Exceptions.h>
#include <DenseVector.h>
#include <DenseMatrix.h>
#include <TwoD_Mapped_Node_Mesh.h>
#include <Utility.h>

namespace CppNoddy {

  template <typename _Type>
  void TwoD_Mapped_Node_Mesh<_Type>::set_nodes_vars(const std::size_t nodex, const std::size_t nodey, const DenseVector<_Type>& U) {
#ifdef PARANOID
    if(U.size() > m_nv) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.set_nodes_vars method is trying to use a \n";
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
  DenseVector<_Type> TwoD_Mapped_Node_Mesh<_Type>::get_nodes_vars(const std::size_t nodex, const std::size_t nodey) const {
#ifdef PARANOID
    if(nodex > m_nx - 1 || nodey > m_ny - 1) {
      std::string problem;
      problem = " The TwoD_Mapped_Node_Mesh.get_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw ExceptionRange(problem, m_nx, nodex, m_ny, nodey);
    }
#endif
    // construct a vector with m_nv elements starting from a pointer
    DenseVector<_Type> nodes_vars(m_nv, &m_vars[(nodex * m_ny + nodey) * m_nv ]);
    return nodes_vars;
  }

  template <typename _Type>
  void TwoD_Mapped_Node_Mesh<_Type>::init_mapping() {
    #ifdef DEBUG
    std::cout << "[DEBUG] Physical domain is [" << m_left << "," << m_right
              << "] x [" << m_bottom << "," << m_top << "]\n";
    std::cout << "[DEBUG] Computational domain is ["
              << FnComp_X(m_left) << "," << FnComp_X(m_right)
              << "] x [" << FnComp_Y(m_bottom) << "," << FnComp_Y(m_top) << "]\n";
    #endif
    {
      // a uniform mesh in the computational coordinates
      double comp_left(FnComp_X(m_left));
      double comp_right(FnComp_X(m_right));
      const double comp_delta_x = (comp_right - comp_left) / (m_nx - 1);
      for(std::size_t i = 0; i < m_nx; ++i) {
        m_compX[i] = comp_left + comp_delta_x * i;
      }
    }
    {
      // a uniform mesh in the computational coordinates
      double comp_bottom(FnComp_Y(m_bottom));
      double comp_top(FnComp_Y(m_top));
      const double comp_delta_y = (comp_top - comp_bottom) / (m_ny - 1);
      for(std::size_t j = 0; j < m_ny; ++j) {
        m_compY[j] = comp_bottom + comp_delta_y * j;
      }
    }

    // temporary fill to span the domain -- corrected below
    m_X = Utility::uniform_node_vector(m_left,m_right,m_nx);
    m_Y = Utility::uniform_node_vector(m_bottom,m_top,m_ny);

    // we now need the corresponding m_Y coordinates in the physical domain
    {
      for(unsigned j = 1; j < m_ny; ++j) {
        // for each node in m_compY
        unsigned kmin(0);
        double min(99e9);
        for(unsigned k = 0; k < m_ny; ++k) {
          // find the y value that is closest to it
          if(std::abs(FnComp_Y(m_Y[k]) - m_compY[j]) < min) {
            min = std::abs(FnComp_Y(m_Y[k]) - m_compY[j]);
            kmin = k;
          }
        }
        double y = m_Y[kmin];
        double delta = 1.e-8;
        double correction = 1.0;
        do {
          double newY = FnComp_Y(y + delta) - m_compY[j];
          double oldY = FnComp_Y(y) - m_compY[j];
          double deriv = (newY-oldY)/delta;
          correction = -oldY/deriv;
          y += correction;
        } while(fabs(correction) > 1.e-8);
        m_Y[ j ] = y;
      }
    }
    // we now need the corresponding m_X coordinates in the physical domain
    {
      for(unsigned i = 1; i < m_nx; ++i) {
        // for each node in m_compY
        unsigned kmin(0);
        double min(99e9);
        for(unsigned k = 0; k < m_ny; ++k) {
          // find the y value that is closest to it
          if(std::abs(FnComp_X(m_X[k]) - m_compX[i]) < min) {
            min = std::abs(FnComp_X(m_X[k]) - m_compX[i]);
            kmin = k;
          }
        }
        double x = m_X[kmin];
        double delta = 1.e-8;
        double correction = 1.0;
        do {
          double newX = FnComp_X(x + delta) - m_compX[i];
          double oldX = FnComp_X(x) - m_compX[i];
          double deriv = (newX-oldX)/delta;
          correction = -oldX/deriv;
          x += correction;
        } while(fabs(correction) > 1.e-8);
        m_X[ i ] = x;
      }
    }
  }

  template <typename _Type>
  std::pair< double, double> TwoD_Mapped_Node_Mesh<_Type>::get_comp_step_sizes() const {
    std::pair<double,double> steps;
    steps.first = m_compX[1]-m_compX[0];
    steps.second = m_compY[1]-m_compY[0];
    return steps;
  }

  template <typename _Type>
  std::pair< std::size_t, std::size_t > TwoD_Mapped_Node_Mesh<_Type>::get_nnodes() const {
    std::pair< std::size_t, std::size_t > nodes;
    nodes.first = m_nx;
    nodes.second = m_ny;
    return nodes;
  }

  template <typename _Type>
  std::size_t TwoD_Mapped_Node_Mesh<_Type>::get_nvars() const {
    return m_nv;
  }

  template <typename _Type>
  DenseVector<double>& TwoD_Mapped_Node_Mesh<_Type>::xnodes() {
    return m_X;
  }

  template <typename _Type>
  DenseVector<double>& TwoD_Mapped_Node_Mesh<_Type>::ynodes() {
    return m_Y;
  }

  template<>
  void TwoD_Mapped_Node_Mesh<double>::normalise(const std::size_t& var) {
    double maxval(max(var));
    m_vars.scale(1./maxval);
  }

  template<>
  void TwoD_Mapped_Node_Mesh<D_complex>::normalise(const std::size_t& var) {
    //std::cout << "[DEBUG] asked to normalise a complex mesh\n";
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
    //std::cout << "[DEBUG] MAX |variable| had complex value of " << factor << "\n";
    m_vars.scale(1./factor);
  }

  template<>
  void TwoD_Mapped_Node_Mesh<double>::dump_gnu(std::string filename) const {
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
  void TwoD_Mapped_Node_Mesh<D_complex>::dump_gnu(std::string filename) const {
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

  //the templated versions we require are:
  template class TwoD_Mapped_Node_Mesh<double>
  ;
  template class TwoD_Mapped_Node_Mesh<std::complex<double> >
  ;

}
