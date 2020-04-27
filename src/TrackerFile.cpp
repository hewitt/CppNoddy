/// \file TrackerFile.cpp
/// Implementation details for the TrackerFile object.

#include <iomanip>

#include <TrackerFile.h>
#include <Utility.h>
#include <OneD_Node_Mesh.h>

namespace CppNoddy {

  TrackerFile::TrackerFile(int prec) {
    this -> m_precision = prec;
  }

  TrackerFile::TrackerFile(std::string filename, int prec) {
    m_dumpFile.open(filename.c_str());
    this -> m_precision = prec;
    m_dumpFile.precision(prec);
    m_dumpFile << std::scientific;
  }

  TrackerFile::~TrackerFile() {
    m_dumpFile.close();
  }

  void TrackerFile::push_ptr(double* scalar, std::string desc) {
    p_doubles.push_back(scalar);
    m_doubleDesc.push_back(desc);
  }

  void TrackerFile::push_ptr(D_complex* scalar, std::string desc) {
    //p_Doubles.push_back( &( scalar -> real() ) );
    //p_Doubles.push_back( &( scalar -> imag() ) );
    // above is no longer valid for g++-7.2
    p_doubles.push_back(&reinterpret_cast<double(&)[2]>(*scalar)[0]);
    p_doubles.push_back(&reinterpret_cast<double(&)[2]>(*scalar)[1]);
    m_doubleDesc.push_back(desc + " (real)");
    m_doubleDesc.push_back(desc + " (imag)");
  }

  void TrackerFile::push_ptr(DenseVector<double>* ptr_to_vector, std::string desc) {
    p_doubleVectors.push_back(ptr_to_vector);
    m_doubleVectorDesc.push_back(desc);
  }

  void TrackerFile::push_ptr(DenseVector<D_complex>* ptr_to_vector, std::string desc) {
    p_complexVectors.push_back(ptr_to_vector);
    m_complexVectorDesc.push_back(desc);
  }

  void TrackerFile::push_ptr(OneD_Node_Mesh<double>* ptr_to_mesh, std::string desc) {
    p_doubleMesh.push_back(ptr_to_mesh);
    m_doubleMeshDesc.push_back(desc);
  }

  void TrackerFile::push_ptr(OneD_Node_Mesh<D_complex>* ptr_to_mesh, std::string desc) {
    p_complexMesh.push_back(ptr_to_mesh);
    m_complexMeshDesc.push_back(desc);
  }

  void TrackerFile::push_ptr(OneD_Node_Mesh<D_complex, D_complex>* ptr_to_mesh, std::string desc) {
    p_complexComplexMesh.push_back(ptr_to_mesh);
    m_complexMeshDesc.push_back(desc);
  }
  void TrackerFile::newline() {
    m_dumpFile << "\n";
  }

  void TrackerFile::set_filename(std::string filename) {
    m_dumpFile.close();
    m_dumpFile.open(filename.c_str());
    m_dumpFile.precision(m_precision);
  }

  void TrackerFile::precision(unsigned prec) {
    this -> m_precision = prec;
    m_dumpFile.precision(prec);
  }

  void TrackerFile::header() {
    // write the header
    m_dumpFile << " # Header : \n # ";
    for(std::size_t i = 0; i < m_doubleDesc.size(); ++i) {
      m_dumpFile << m_doubleDesc[ i ] << " | ";
    }
    for(std::size_t i = 0; i < m_doubleVectorDesc.size(); ++i) {
      m_dumpFile << m_doubleVectorDesc[ i ] << " | ";
    }
    for(std::size_t i = 0; i < m_complexVectorDesc.size(); ++i) {
      m_dumpFile << m_complexVectorDesc[ i ] + " (Real)" << " | ";
      m_dumpFile << m_complexVectorDesc[ i ] + " (Imag)" << " | ";
    }
    for(std::size_t i = 0; i < m_doubleMeshDesc.size(); ++i) {
      m_dumpFile << m_doubleMeshDesc[ i ] + " (nodes)" << " | ";
      for(unsigned var = 0; var < (*p_doubleMesh[ i ]).get_nvars(); ++var) {
        m_dumpFile << m_doubleMeshDesc[ i ] + " var#" + Utility::stringify(var) + " | ";
      }
    }
    for(std::size_t i = 0; i < m_complexMeshDesc.size(); ++i) {
      m_dumpFile << m_complexMeshDesc[ i ] + "(nodes)" << " | ";
      for(unsigned var = 0; var < (*p_complexMesh[ i ]).get_nvars(); ++var) {
        m_dumpFile << m_complexMeshDesc[ i ] + " var#" + Utility::stringify(var) + " (Real) | ";
        m_dumpFile << m_complexMeshDesc[ i ] + " var#" + Utility::stringify(var) + " (Imag) | ";
      }
    }
    for(std::size_t i = 0; i < m_complexComplexMeshDesc.size(); ++i) {
      m_dumpFile << m_complexComplexMeshDesc[ i ] + "(nodes)_real" << " | ";
      m_dumpFile << m_complexComplexMeshDesc[ i ] + "(nodes)_imag" << " | ";
      for(unsigned var = 0; var < (*p_complexComplexMesh[ i ]).get_nvars(); ++var) {
        m_dumpFile << m_complexComplexMeshDesc[ i ] + " var#" + Utility::stringify(var) + " (Real) | ";
        m_dumpFile << m_complexComplexMeshDesc[ i ] + " var#" + Utility::stringify(var) + " (Imag) | ";
      }
    }
    m_dumpFile << "\n";
  }

  void TrackerFile::update() {
    unsigned block_size(1);
    if(!p_doubleVectors.empty()) {
      block_size = p_doubleVectors[ 0 ] -> size();
    }
    if(!p_complexVectors.empty()) {
      block_size = p_complexVectors[ 0 ] -> size();
    }
    if(!p_doubleMesh.empty()) {
      block_size = p_doubleMesh[ 0 ] -> get_nnodes();
    }
    if(!p_complexMesh.empty()) {
      block_size = p_complexMesh[ 0 ] -> get_nnodes();
    }
    if(!p_complexComplexMesh.empty()) {
      block_size = p_complexComplexMesh[ 0 ] -> get_nnodes();
    }

    for(unsigned line = 0; line < block_size; ++line) {
      dump_scalar_data();
      if(!p_doubleVectors.empty()) {
        // for each vector ptr
        for(std::size_t i = 0; i < p_doubleVectors.size(); ++i) {
          m_dumpFile << (*p_doubleVectors[ i ]) [ line ] << " ";
        }
      }
      if(!p_complexVectors.empty()) {
        // for each vector ptr
        for(std::size_t i = 0; i < p_complexVectors.size(); ++i) {
          m_dumpFile << (*p_complexVectors[ i ]) [ line ].real() << " ";
          m_dumpFile << (*p_complexVectors[ i ]) [ line ].imag() << " ";
        }
      }
      if(!p_doubleMesh.empty()) {
        // for each mesh ptr
        for(std::size_t i = 0; i < p_doubleMesh.size(); ++i) {
          m_dumpFile << (*p_doubleMesh[ i ]).coord(line) << " ";
          for(unsigned var = 0; var < p_doubleMesh[ i ] -> get_nvars(); ++var) {
            m_dumpFile << (*p_doubleMesh[ i ])(line, var) << " ";
          }
        }
      }
      if(!p_complexMesh.empty()) {
        // for each mesh ptr
        for(std::size_t i = 0; i < p_complexMesh.size(); ++i) {
          m_dumpFile << (*p_complexMesh[ i ]).coord(line) << " ";
          for(unsigned var = 0; var < p_complexMesh[ i ] -> get_nvars(); ++var) {
            m_dumpFile << (*p_complexMesh[ i ])(line, var).real() << " ";
            m_dumpFile << (*p_complexMesh[ i ])(line, var).imag() << " ";
          }
        }
      }
      if(!p_complexComplexMesh.empty()) {
        // for each mesh ptr
        for(std::size_t i = 0; i < p_complexComplexMesh.size(); ++i) {
          m_dumpFile << (*p_complexComplexMesh[ i ]).coord(line).real() << " ";
          m_dumpFile << (*p_complexComplexMesh[ i ]).coord(line).imag() << " ";
          for(unsigned var = 0; var < p_complexComplexMesh[ i ] -> get_nvars(); ++var) {
            m_dumpFile << (*p_complexComplexMesh[ i ])(line, var).real() << " ";
            m_dumpFile << (*p_complexComplexMesh[ i ])(line, var).imag() << " ";
          }
        }
      }
      m_dumpFile << "\n";
    }
    // flush the buffer
    m_dumpFile.flush();
  }

  void TrackerFile::dump_scalar_data() {
    if(!p_doubles.empty()) {
      // simple flat data file
      for(std::size_t i = 0; i < p_doubles.size(); ++i) {
        m_dumpFile << *p_doubles[ i ] << " ";
      }
    }
  }

} //end namepsace
