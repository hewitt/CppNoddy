/// \file DistributedLinearSystem.h
/// Specification of a sparse-storage distributed linear system class.

#if defined(PETSC_Z) || defined(PETSC_D)

#ifndef DISTRIBUTEDLINEARSYSTEM_H
#define DISTRIBUTEDLINEARSYSTEM_H

#include <DistributedMatrix.h>
#include <DistributedVector.h>
#include <Exceptions.h>
#include <vector>

namespace CppNoddy {

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for distributed problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class DistributedLinearSystem {

   public:

    /// Constructor for a distributed linear system object.
    /// \param pA A pointer to the 'A matrix', an NxN double/complex distributed matrix
    /// \param pB A pointer to the 'B vector' a size N double/complex distributed vector
    DistributedLinearSystem(DistributedMatrix<_Type>* pA, DistributedVector<_Type>* pB  ) {
      m_pA = pA;
      m_pB = pB;
      KSPCreate(PETSC_COMM_WORLD,&m_ksp);
      // default the relative tolerance "rtol" to 1.e-10 for iterative solvers
      // KSPSetTolerances( m_ksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    }

    /// Destructor for a linear system object.
    ~DistributedLinearSystem(){
      //VecDestroy(&m_x);
      KSPDestroy(&m_ksp);
    }

    /// Solve the sparse system in place with soln overwriting B
    void solve() {
      // construct using the PETSc matrix
      KSPSetOperators(m_ksp,*(m_pA->get_pMat()),*(m_pA->get_pMat()));

      /////////////////////////////
      // default solver is MUMPS //
      /////////////////////////////
      /*
      KSPSetType(m_ksp,KSPPREONLY);
      // preconditioner
      KSPGetPC(m_ksp,&m_pc);
      // hardwire a DIRECT SOLVER via MUMPS
      PCSetType(m_pc,PCLU);
      PCFactorSetMatSolverType(m_pc,MATSOLVERMUMPS);
      PCFactorSetUpMatSolverType(m_pc);
      */      
      /////////////////////////////
      
      
      // override the KSP member data using command line options
      KSPSetFromOptions(m_ksp);
      KSPSetUp(m_ksp);
      // We'll hard wire the solution back into the RHS vector
      KSPSolve(m_ksp,*(m_pB->get_pVec()),*(m_pB->get_pVec()));
      KSPView(m_ksp, PETSC_VIEWER_STDOUT_WORLD);
    }

    void view() {
      KSPView(m_ksp, PETSC_VIEWER_STDOUT_WORLD);
    }
    
    /// Get a sequential vector from the distributed RHS vector
    /// \return A sequential PETSc Vec object -- obviously not deleted
    /// when the DistributedLinearSystem object destructs.
    Vec get_seq_solution() {
      // context for the scatter
      VecScatter ctx;
      Vec v_seq;
      // manual says no need for VecCreate
      VecScatterCreateToAll(*(m_pB->get_pVec()), &ctx, &v_seq);
      // scatter the data to v_seq
      VecScatterBegin(ctx, *(m_pB->get_pVec()), v_seq, INSERT_VALUES, SCATTER_FORWARD);
      VecScatterEnd(ctx, *(m_pB->get_pVec()), v_seq, INSERT_VALUES, SCATTER_FORWARD);
      // clean up
      VecScatterDestroy(&ctx);
      return v_seq;
    }

    /// Get the solution as a CppNoddy::DenseVector<_Type> vector
    /// \return The entire solution vector as a sequential DenseVector<_Type>
    DenseVector<_Type> get_solution() {
      // convert the distributed vec solution to a sequential one
      Vec seq_soln = get_seq_solution();
      // this "array" is getting a pointer not copying data
      PetscScalar* array;
      VecGetArray(seq_soln,&array);
      // now copy to the CppNoddy densevctor using the constructor
      PetscInt size(0);
      VecGetSize(seq_soln,&size);
      DenseVector<_Type> soln(size,array);
      return soln;
    }
    
    /*
   /// Factorise the Ax=B system
   void factorise(){
   // construct using the PETSc matrix
   KSPSetOperators(m_ksp,*(m_pA->get_pMat()),*(m_pA->get_pMat()));
   //
   KSPSetType(m_ksp,KSPPREONLY);
   // preconditioner
   KSPGetPC(m_ksp,&m_pc);
   // hardwire a DIRECT SOLVER via MUMPS
   PCSetType(m_pc,PCLU);
   PCFactorSetMatSolverType(m_pc,MATSOLVERMUMPS);
   PCFactorSetUpMatSolverType(m_pc);
   //call MatGetFactor() to create F
   PCFactorGetMatrix(m_pc,&m_F);
   // allow override using command line options?
   KSPSetFromOptions(m_ksp);
   KSPSetUp(m_ksp);
   }
   
   /// Resolve the same system using the same factorisation
   void solve_using_factorisation();
    */
    
   private:
    // pointer to a distributed LHS matrix
    DistributedMatrix<_Type>* m_pA;
    // pointer to the distributed RHS vector
    DistributedVector<_Type>* m_pB;
    // Solver
    KSP m_ksp;
    // Preconditioner
    PC m_pc;
    // factors
    Mat m_F;
  };

} //end namepsace

#endif // include guard

#endif // PETSc_D or PETSc_Z
