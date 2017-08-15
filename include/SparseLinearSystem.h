/// \file SparseLinearSystem.h
/// Specification of a sparse-storage linear system class.

#ifndef SPARSELINEARSYSTEM_H
#define SPARSELINEARSYSTEM_H

#include <SparseMatrix.h>
#include <DenseVector.h>
#include <Exceptions.h>
#include <LinearSystem_base.h>


#ifdef SUPERLU
  // There are conflicts between Superlu3 include files
  // slu_ddefs.h & slu_zdefs.h for double and complex systems.
  // Hence SLU.h is these two includes put together, but then
  // separated by the SLUD and SLUZ namespaces.
  #include <SLU.h>
#endif

#ifdef MUMPS_SEQ
  // double precision real includes
  #include "dmumps_c.h"
  // double precision complex includes
  #include "zmumps_c.h"
#endif

namespace CppNoddy
{
  /// A linear system class for vector right-hand sides.
  /// The class is constructed for SPARSE problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class SparseLinearSystem : public LinearSystem_base
  {

  public:

    /// Constructor for a sparse linear system object.
    /// \param Aptr A pointer to the 'A matrix', an NxN double/complex sparse matrix
    /// \param Bptr A pointer to the 'B vector' a size N double/complex dense vector
    /// \param which A string that indicates which solver to use: native (default), superlu or mumps_seq 
    SparseLinearSystem( SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which = "native" );

    /// Destructor for a linear system object.
    ~SparseLinearSystem()
    {}

    /// Solve the sparse system
    void solve();

  private:

    /// Solve the linear system using the native elimination -- quite a naive implementation I imagine.
    void solve_native();

    /// Solve the linear system by linking to the SuperLU library
    void solve_superlu();

    /// Solve the linear system by linking to the MUMPS (Sequential) library
    void solve_mumps_seq();

    /// Back substitution routine for dense systems.
    /// \param A The upper triangular matrix LHS
    /// \param B The dense vector RHS
    void backsub( SparseMatrix<_Type> &A, DenseVector<_Type> &B ) const;

    /// check on pivot size for the native elimination routine
    double MIN_PIV;

    /// pointer to a sparse LHS matrix
    SparseMatrix<_Type>* p_A;
    /// pointer to the RHS vector
    DenseVector<_Type>* p_B;

  };

} //end namepsace
#endif
