/// \file TrackerFile.cpp
/// Implementation details for the TrackerFile object.

#include <iomanip>

#include <TrackerFile.h>
#include <Utility.h>
#include <OneD_Node_Mesh.h>

namespace CppNoddy {

  TrackerFile::TrackerFile(int prec) {
    this -> PREC = prec;
  }

  TrackerFile::TrackerFile(std::string filename, int prec) {
    dumpfile.open(filename.c_str());
    this -> PREC = prec;
    dumpfile.precision(prec);
    dumpfile << std::scientific;
  }

  TrackerFile::~TrackerFile() {
    dumpfile.close();
  }

  void TrackerFile::push_ptr(double* scalar, std::string desc) {
    p_DOUBLES.push_back(scalar);
    DOUBLE_DESC.push_back(desc);
  }

  void TrackerFile::push_ptr(D_complex* scalar, std::string desc) {
    //p_DOUBLES.push_back( &( scalar -> real() ) );
    //p_DOUBLES.push_back( &( scalar -> imag() ) );
    // above is no longer valid for g++-7.2
    p_DOUBLES.push_back(&reinterpret_cast<double(&)[2]>(*scalar)[0]);
    p_DOUBLES.push_back(&reinterpret_cast<double(&)[2]>(*scalar)[1]);
    DOUBLE_DESC.push_back(desc + " (real)");
    DOUBLE_DESC.push_back(desc + " (imag)");
  }

  void TrackerFile::push_ptr(DenseVector<double>* ptr_to_vector, std::string desc) {
    p_DVECTORS.push_back(ptr_to_vector);
    DVECTOR_DESC.push_back(desc);
  }

  void TrackerFile::push_ptr(DenseVector<D_complex>* ptr_to_vector, std::string desc) {
    p_CVECTORS.push_back(ptr_to_vector);
    CVECTOR_DESC.push_back(desc);
  }

  void TrackerFile::push_ptr(OneD_Node_Mesh<double>* ptr_to_mesh, std::string desc) {
    p_DMESH.push_back(ptr_to_mesh);
    DMESH_DESC.push_back(desc);
  }

  void TrackerFile::push_ptr(OneD_Node_Mesh<D_complex>* ptr_to_mesh, std::string desc) {
    p_CMESH.push_back(ptr_to_mesh);
    CMESH_DESC.push_back(desc);
  }

  void TrackerFile::push_ptr(OneD_Node_Mesh<D_complex, D_complex>* ptr_to_mesh, std::string desc) {
    p_CCMESH.push_back(ptr_to_mesh);
    CCMESH_DESC.push_back(desc);
  }
  void TrackerFile::newline() {
    dumpfile << "\n";
  }

  void TrackerFile::set_filename(std::string filename) {
    dumpfile.close();
    dumpfile.open(filename.c_str());
    dumpfile.precision(PREC);
  }

  void TrackerFile::precision(unsigned prec) {
    this -> PREC = prec;
    dumpfile.precision(prec);
  }

  void TrackerFile::header() {
    // write the header
    dumpfile << " # Header : \n # ";
    for(std::size_t i = 0; i < DOUBLE_DESC.size(); ++i) {
      dumpfile << DOUBLE_DESC[ i ] << " | ";
    }
    for(std::size_t i = 0; i < DVECTOR_DESC.size(); ++i) {
      dumpfile << DVECTOR_DESC[ i ] << " | ";
    }
    for(std::size_t i = 0; i < CVECTOR_DESC.size(); ++i) {
      dumpfile << CVECTOR_DESC[ i ] + " (Real)" << " | ";
      dumpfile << CVECTOR_DESC[ i ] + " (Imag)" << " | ";
    }
    for(std::size_t i = 0; i < DMESH_DESC.size(); ++i) {
      dumpfile << DMESH_DESC[ i ] + " (nodes)" << " | ";
      for(unsigned var = 0; var < (*p_DMESH[ i ]).get_nvars(); ++var) {
        dumpfile << DMESH_DESC[ i ] + " var#" + Utility::stringify(var) + " | ";
      }
    }
    for(std::size_t i = 0; i < CMESH_DESC.size(); ++i) {
      dumpfile << CMESH_DESC[ i ] + "(nodes)" << " | ";
      for(unsigned var = 0; var < (*p_CMESH[ i ]).get_nvars(); ++var) {
        dumpfile << CMESH_DESC[ i ] + " var#" + Utility::stringify(var) + " (Real) | ";
        dumpfile << CMESH_DESC[ i ] + " var#" + Utility::stringify(var) + " (Imag) | ";
      }
    }
    for(std::size_t i = 0; i < CCMESH_DESC.size(); ++i) {
      dumpfile << CCMESH_DESC[ i ] + "(nodes)_real" << " | ";
      dumpfile << CCMESH_DESC[ i ] + "(nodes)_imag" << " | ";
      for(unsigned var = 0; var < (*p_CCMESH[ i ]).get_nvars(); ++var) {
        dumpfile << CCMESH_DESC[ i ] + " var#" + Utility::stringify(var) + " (Real) | ";
        dumpfile << CCMESH_DESC[ i ] + " var#" + Utility::stringify(var) + " (Imag) | ";
      }
    }
    dumpfile << "\n";
  }

  void TrackerFile::update() {
    unsigned block_size(1);
    if(!p_DVECTORS.empty()) {
      block_size = p_DVECTORS[ 0 ] -> size();
    }
    if(!p_CVECTORS.empty()) {
      block_size = p_CVECTORS[ 0 ] -> size();
    }
    if(!p_DMESH.empty()) {
      block_size = p_DMESH[ 0 ] -> get_nnodes();
    }
    if(!p_CMESH.empty()) {
      block_size = p_CMESH[ 0 ] -> get_nnodes();
    }
    if(!p_CCMESH.empty()) {
      block_size = p_CCMESH[ 0 ] -> get_nnodes();
    }

    for(unsigned line = 0; line < block_size; ++line) {
      dump_scalar_data();
      if(!p_DVECTORS.empty()) {
        // for each vector ptr
        for(std::size_t i = 0; i < p_DVECTORS.size(); ++i) {
          dumpfile << (*p_DVECTORS[ i ]) [ line ] << " ";
        }
      }
      if(!p_CVECTORS.empty()) {
        // for each vector ptr
        for(std::size_t i = 0; i < p_CVECTORS.size(); ++i) {
          dumpfile << (*p_CVECTORS[ i ]) [ line ].real() << " ";
          dumpfile << (*p_CVECTORS[ i ]) [ line ].imag() << " ";
        }
      }
      if(!p_DMESH.empty()) {
        // for each mesh ptr
        for(std::size_t i = 0; i < p_DMESH.size(); ++i) {
          dumpfile << (*p_DMESH[ i ]).coord(line) << " ";
          for(unsigned var = 0; var < p_DMESH[ i ] -> get_nvars(); ++var) {
            dumpfile << (*p_DMESH[ i ])(line, var) << " ";
          }
        }
      }
      if(!p_CMESH.empty()) {
        // for each mesh ptr
        for(std::size_t i = 0; i < p_CMESH.size(); ++i) {
          dumpfile << (*p_CMESH[ i ]).coord(line) << " ";
          for(unsigned var = 0; var < p_CMESH[ i ] -> get_nvars(); ++var) {
            dumpfile << (*p_CMESH[ i ])(line, var).real() << " ";
            dumpfile << (*p_CMESH[ i ])(line, var).imag() << " ";
          }
        }
      }
      if(!p_CCMESH.empty()) {
        // for each mesh ptr
        for(std::size_t i = 0; i < p_CCMESH.size(); ++i) {
          dumpfile << (*p_CCMESH[ i ]).coord(line).real() << " ";
          dumpfile << (*p_CCMESH[ i ]).coord(line).imag() << " ";
          for(unsigned var = 0; var < p_CCMESH[ i ] -> get_nvars(); ++var) {
            dumpfile << (*p_CCMESH[ i ])(line, var).real() << " ";
            dumpfile << (*p_CCMESH[ i ])(line, var).imag() << " ";
          }
        }
      }
      dumpfile << "\n";
    }
    // flush the buffer
    dumpfile.flush();
  }

  void TrackerFile::dump_scalar_data() {
    if(!p_DOUBLES.empty()) {
      // simple flat data file
      for(std::size_t i = 0; i < p_DOUBLES.size(); ++i) {
        dumpfile << *p_DOUBLES[ i ] << " ";
      }
    }
  }

} //end namepsace
