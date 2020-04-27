/// \file TrackerFile.h
/// A class that can be passed pointers to scalar/vector/mesh
/// objects, then their contents are written to a file each
/// time that the update method is called.

#ifndef TRACKERFILE_H
#define TRACKERFILE_H

#include <fstream>
#include <iostream>
#include <string>

#include <Types.h>
#include <OneD_Node_Mesh.h>

namespace CppNoddy {
  class TrackerFile {
   public:

    TrackerFile(int prec = 12);

    TrackerFile(std::string filename, int prec = 12);

    ~TrackerFile();

    void precision(unsigned prec);

    void push_ptr(double* scalar, std::string desc = "");

    void push_ptr(D_complex* scalar, std::string desc = "");

    void push_ptr(DenseVector<double>* ptr_to_vector, std::string desc = "");

    void push_ptr(DenseVector<D_complex>* ptr_to_vector, std::string desc = "");

    void push_ptr(OneD_Node_Mesh<double>* ptr_to_mesh, std::string desc = "");

    void push_ptr(OneD_Node_Mesh<D_complex>* ptr_to_mesh, std::string desc = "");

    void push_ptr(OneD_Node_Mesh<D_complex, D_complex>* ptr_to_mesh, std::string desc = "");

    void newline();

    void set_filename(std::string filename);

    void header();

    void update();

   protected:

    mutable std::ofstream m_dumpFile;

   private:

    /// private output method to write all scalar data to dumpfile
    void dump_scalar_data();

    /// a vector of description strings
    std::vector< std::string > m_doubleDesc;
    std::vector< std::string > m_doubleVectorDesc;
    std::vector< std::string > m_complexVectorDesc;
    std::vector< std::string > m_doubleMeshDesc;
    std::vector< std::string > m_complexMeshDesc;
    std::vector< std::string > m_complexComplexMeshDesc;
    /// a vector of pointers to scalars
    std::vector< double* > p_doubles;
    /// a vector of pointers to double dense vectors
    std::vector< DenseVector<double>* > p_doubleVectors;
    /// a vector of pointers to complex dense vectors
    std::vector< DenseVector<D_complex>* > p_complexVectors;
    /// a vector of pointers to real 1D meshes
    std::vector< OneD_Node_Mesh<double>* > p_doubleMesh;
    /// a vector of pointers to complex 1D meshes
    std::vector< OneD_Node_Mesh<D_complex>* > p_complexMesh;
    /// a vector of pointers to complex-complex 1D meshes
    std::vector< OneD_Node_Mesh<D_complex, D_complex>* > p_complexComplexMesh;

    /// output precision
    unsigned m_precision;

  }
  ; //end class
} //end namepsace
#endif
