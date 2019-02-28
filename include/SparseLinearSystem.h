/// \file SparseLinearSystem.h
/// Specification of a sparse-storage linear system class.

#ifndef SPARSELINEARSYSTEM_H
#define SPARSELINEARSYSTEM_H

#include <SparseMatrix.h>
#include <DenseVector.h>
#include <Exceptions.h>

#if defined(PETSC_Z) || defined(PETSC_D)
#include "petscksp.h"
#endif

namespace CppNoddy {
  /// A linear system class for vector right-hand sides.
  /// The class is constructed for SPARSE problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class SparseLinearSystem {

   public:

    /// Constructor for a sparse linear system object.
    /// \param Aptr A pointer to the 'A matrix', an NxN double/complex sparse matrix
    /// \param Bptr A pointer to the 'B vector' a size N double/complex dense vector
    /// \param which A string that indicates which solver to use: native (default) pr petsc
    SparseLinearSystem(SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which = "native");

    /// Destructor for a linear system object.
    ~SparseLinearSystem();

    /// deallocates some objects
    void cleanup();

    /// Solve the sparse system
    void solve();

    /// Factorise the Ax=B system
    void factorise();

    /// Resolve the same system using the same factorisation
    void solve_using_factorisation();

    void temp_solve();

   private:
    // solve by linking to PETSc
    void solve_petsc();

    // factorise by linking to PETSc -- allow for re-solves
    void factorise_petsc();
    
    /// pointer to a sparse LHS matrix
    SparseMatrix<_Type>* m_pA;
    /// pointer to the RHS vector
    DenseVector<_Type>* m_pB;

    /// a string ID to pick out the appropriate solver
    std::string m_version;

    /// indicates that the matrix has been factorised
    bool m_factorised;

#if defined(PETSC_Z) || defined(PETSC_D)
    Vec            m_petsc_x,m_petsc_B;   /* B = RHS and x = soln */
    Mat            m_petsc_F;
    KSP            m_petsc_ksp;     /* linear solver context */
    PC             m_petsc_pc;      /* preconditioner -- though hard wired for MUMPS direct method */
    PetscMPIInt    m_petsc_rank, m_petsc_size;
#endif
  };

} //end namepsace
#endif
