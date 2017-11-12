/// \file SparseLinearSystem.h
/// Specification of a sparse-storage linear system class.

#ifndef SPARSELINEARSYSTEM_H
#define SPARSELINEARSYSTEM_H

#include <SparseMatrix.h>
#include <DenseVector.h>
#include <Exceptions.h>
#include <LinearSystem_base.h>

#if defined(PETSC_Z) || defined(PETSC_D)
  #include <petscksp.h>
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
    /// \param which A string that indicates which solver to use: native (default) pr petsc
    SparseLinearSystem( SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which = "native" );

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

  private:
    // solve by linking to PETSc
    void solve_petsc();

    // factorise by linking to PETSc -- allow for re-solves
    void factorise_petsc();

    /// Solve the linear system using the native elimination -- quite a naive implementation I imagine.
    void solve_native();

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

    /// indicates that the matrix has been factorised
    bool factorised_;

    #if defined(PETSC_Z) || defined(PETSC_D)
      Vec            x_,B_;      /* B = RHS and x = soln */
      Mat            F_;
      KSP            ksp_;       /* linear solver context */
      PC             pc_;        /* preconditioner -- though hard wired for MUMPS direct method */
      PetscMPIInt    rank_, size_;
    #endif
  };

} //end namepsace
#endif
