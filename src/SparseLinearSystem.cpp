/// \file SparseLinearSystem.cpp
/// Implementation for the LinearSystem class

#include <vector>
#include <set>

#include <SparseLinearSystem.h>
#include <Exceptions.h>
#include <Types.h>
#include <Timer.h>


#if defined(PETSC_D) || defined(PETSC_Z)
#include "petsc.h"
#include "mpi.h"
#endif

namespace CppNoddy {

  template <typename _Type>
  SparseLinearSystem<_Type>::SparseLinearSystem(SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which) : m_factorised(false) {
    m_pA = Aptr;
    m_pB = Bptr;
    m_version = which;
    //
    if (m_version != "petsc") {
      std::string problem;
      problem = "The SparseLinearSystem has been instantiated with an unrecognised\n";
      problem += "request for a solver type. Options: 'petsc'. \n";
      throw ExceptionRuntime(problem);
    }
    //
    if(m_version == "petsc") {
#if !defined(PETSC_D) && !defined(PETSC_Z)
      std::string problem;
      problem = "The SparseLinearSystem has been instantiated for a petsc solver.\n";
      problem += "HOWEVER, CppNoddy was not compiled with PETSC_D or PETSC_Z defined.\n";
      throw ExceptionRuntime(problem);
#endif
    }
    //
#if defined(PETSC_D) || defined(PETSC_Z)
    if(m_version == "petsc") {
      int flag(0);
      MPI_Initialized(&flag);
      if(flag != 1) {
        std::string problem;
        problem = "The SparseLinearSystem has been instantiated for a petsc solver.\n";
        problem += "You must call PetscInitialize before calling the petsc solver.\n";
        throw ExceptionRuntime(problem);
      }
    }
    MPI_Comm_size(MPI_COMM_WORLD,&m_petsc_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&m_petsc_rank);
    if(m_petsc_size > 1) {
      std::string problem;
      problem = " The SparseLinearSystem/SparseMatrix objects are not written\n";
      problem += " to make use of MPI. Try the DistributedLinearsystem and\n";
      problem += " DistributedMatrix objects if MPI is to be used.\n";
      throw ExceptionRuntime(problem);
    }
#endif
  }

  template<typename _Type>
  SparseLinearSystem<_Type>::~SparseLinearSystem() {
    cleanup();
  }

  template<typename _Type>
  void SparseLinearSystem<_Type>::cleanup() {
#if defined(PETSC_D) || defined(PETSC_Z)
    // delete objects used in the factorisation?
    if(m_factorised) {
      VecDestroy(&m_petsc_x);
      VecDestroy(&m_petsc_B);
      KSPDestroy(&m_petsc_ksp);
      m_factorised = false;
    }
#endif
  }

  template <typename _Type>
  void SparseLinearSystem<_Type>::solve() {
    if("petsc" == m_version) {
      std::cout << "[DEBUG] Pre-factorise\n";
      factorise();
      std::cout << "[DEBUG] Post-factorise\n";
      std::cout << "[DEBUG] Pre-solve using factors\n";
      solve_using_factorisation();
      std::cout << "[DEBUG] Post-solve using factors\n";
    } else { // we catch incorrect m_version choices in the ctor
      std::string problem;
      problem = "CppNoddy needs to be linked to PETSc to solve sparse\n";
      problem += "linear systems.";
      throw ExceptionExternal(problem);
    }
  }


  /*
    $PETSC_ARCH points to a DOUBLE implementation routines are below
  */

  template<>
  void SparseLinearSystem<double>::factorise() {
#if !defined(PETSC_D)
    std::string problem;
    problem = "CppNoddy is linked against the COMPLEX version of PETSc\n";
    problem += "but you are trying to factorise a DOUBLE matrix. Either\n";
    problem += "redefine your matrix as complex, or recompile with $PETSC_ARCH\n";
    problem += "pointing to a DOUBLE version of the PETSc code.";
    throw ExceptionExternal(problem);
#endif
#if defined(PETSC_D)

    if(m_factorised) {
      // already factorised -- so delete and re-create below
      cleanup();
    }

    // store a boolean to indicate that we
    m_factorised = true;
    PetscInt Istart,Iend,n;
    Mat petsc_A;

    // size of the (assumed square) matrix
    n = m_pA -> nrows();
    /*
       Create parallel vectors B (for RHS) & x (for soln).
    */
    VecCreate(PETSC_COMM_WORLD,&m_petsc_B);
    VecSetSizes(m_petsc_B,PETSC_DECIDE,m_pA->nrows());
    // add any command line configuration
    VecSetFromOptions(m_petsc_B);
    // make a copy to define x
    VecDuplicate(m_petsc_B,&m_petsc_x);

    
    // configure the A matrix
    MatCreate(PETSC_COMM_WORLD,&petsc_A);
    // set A to be an nxn matrix
    MatSetSizes(petsc_A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    // add any command line configuration
    MatSetFromOptions(petsc_A);

    // get: all_rows_nnz[i] is the number of nonzero elts in row i
    PetscInt* all_rows_nnz = new PetscInt[ n ];
    m_pA -> nelts_all_rows(all_rows_nnz);

    // pre-allocate memory using the number of non-zero elts
    // in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(petsc_A, 0, all_rows_nnz);
    
    // add any command line configuration
    MatSetFromOptions(petsc_A);
    // finish the A definition
    MatSetUp(petsc_A);
    
    /*
       petsc_A is defined as Sequential above, so this is not strictly needed
    */
    MatGetOwnershipRange(petsc_A,&Istart,&Iend);
    // populate the A matrix using the CppNoddy sparse matrix data
    for(PetscInt i = Istart; i<Iend; ++i) {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = all_rows_nnz[i]; 
      // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
      PetscInt* cols = new PetscInt[all_rows_nnz[i]];
      // store the non-zero elts in this row
      PetscScalar* storage = new PetscScalar[all_rows_nnz[i]];
      // get the data from the CppNoddy sparse matrix structure
      m_pA -> get_row_petsc(i, storage, cols);
      MatSetValues(petsc_A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
      // delete temp storage made in the conversion
      delete[] cols;
      delete[] storage;
    }
    
    // MatSetValue inserted values are generally cached
    // so we need to explicitly do final assembly
    MatAssemblyBegin(petsc_A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(petsc_A,MAT_FINAL_ASSEMBLY);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  Create the linear solver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    */
    KSPCreate(PETSC_COMM_WORLD,&m_petsc_ksp);
    KSPSetOperators(m_petsc_ksp,petsc_A,petsc_A);

    
    /////////////////////////////
    // default solver is MUMPS //
    /////////////////////////////
    KSPSetType(m_petsc_ksp,KSPPREONLY);
    KSPGetPC(m_petsc_ksp,&m_petsc_pc);
    // hardwire a DIRECT SOLVER via MUMPS
    PCSetType(m_petsc_pc,PCLU);
    PCFactorSetMatSolverType(m_petsc_pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverType(m_petsc_pc);
    /////////////////////////////
    

    /* create m_petsc_F */
    PCFactorGetMatrix(m_petsc_pc,&m_petsc_F);

    /* increase estimate of working space by 20% */
    // MatMumpsSetIcntl(m_petsc_F,14,20.0);

    /* sequential ordering */
    // icntl = 7;
    // ival = 2;
    // MatMumpsSetIcntl(m_petsc_F,icntl,ival);

    /* threshhold for row pivot detection */
    // MatMumpsSetIcntl(m_petsc_F,24,1);
    // icntl = 3;
    // val = 1.e-6;
    // MatMumpsSetCntl(m_petsc_F,icntl,val);

    /* compute determinant of A */
    // MatMumpsSetIcntl(m_petsc_F,33,1);
    /* not used unless we initialise PETSc using the command line options */
    // KSPSetFromOptions(m_petsc_ksp);
    
    KSPSetFromOptions(m_petsc_ksp);
    KSPSetUp(m_petsc_ksp);

    // /* determinant calculation */
    // {
    //   PetscInt  icntl,infog34;
    //   PetscReal cntl,rinfo12,rinfo13;
    //   icntl = 3;
    //   MatMumpsGetCntl(F,icntl,&cntl);
    //   // output determinant only the first proc.
    //   if (!rank)
    //   {
    //     MatMumpsGetInfog(F,34,&infog34);
    //     MatMumpsGetRinfog(F,12,&rinfo12);
    //     MatMumpsGetRinfog(F,13,&rinfo13);
    //     PetscPrintf(PETSC_COMM_SELF,"  Mumps row pivot threshhold = %g\n",cntl);
    //     PetscPrintf(PETSC_COMM_SELF,"  Mumps determinant = (%g, %g) * 2^%D \n",(double)rinfo12,(double)rinfo13,infog34);
    //   }
    // }

    MatDestroy(&petsc_A);
    delete[] all_rows_nnz;
    
#endif // check for PETSC_D/Z
  }




  template <>
  void SparseLinearSystem<double>::solve_using_factorisation() {
#if !defined(PETSC_D)
    std::string problem;
    problem = "CppNoddy is linked against the COMPLEX version of PETSc\n";
    problem += "but you are trying to solve a DOUBLE matrix. Either\n";
    problem += "redefine your matrix as complex, or recompile with $PETSC_ARCH\n";
    problem += "pointing to a DOUBLE version of the PETSc code.";
    throw ExceptionExternal(problem);
#endif
#if defined(PETSC_D)
    // size of the (assumed square) matrix
    PetscInt n = m_pA -> nrows();

    // populate the RHS vector using the CppNoddy DenseVector content
    for(PetscInt i = 0; i < n; ++i) {
      VecSetValue(m_petsc_B,i,m_pB->operator[](i),INSERT_VALUES);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the linear system
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    KSPSolve(m_petsc_ksp,m_petsc_B,m_petsc_x);

    /* We can now gather the parallel result back to ALL processes
      This is temporary as the SparseMatrix is stored on each processes
      and is too dumb for "proper" parallelization */
    Vec y;
    // a scatter context
    VecScatter ctx = 0;
    // map all elts of the parallel vector to a sequential copy
    VecScatterCreateToAll(m_petsc_x,&ctx,&y);
    // scatter it
    VecScatterBegin(ctx,m_petsc_x,y,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_petsc_x,y,INSERT_VALUES,SCATTER_FORWARD);
    // clean up
    VecScatterDestroy(&ctx);
    // this array is a pointer not a copy
    PetscScalar* array;
    VecGetArray(y,&array);
    // now copy to the CppNoddy densevctor
    for(PetscInt i=0; i<n; i++) {
      m_pB -> operator[](i) = array[i];
    }
    // follow the docs and Restore after get
    VecRestoreArray(m_petsc_x,&array);
    VecDestroy(&y);
#endif
  }






  /*
    $PETSC_ARCH points to a COMPLEX implementation routines are below
  */

  template<>
  void SparseLinearSystem<std::complex<double> >::factorise() {
#if !defined(PETSC_Z)
    std::string problem;
    problem = "CppNoddy is linked against the DOUBLE version of PETSc\n";
    problem += "but you are trying to factorise a COMPLEX matrix.\n";
    problem += "Recompile with $PETSC_ARCH\n";
    problem += "pointing to a COMPLEX version of the PETSc code.";
    throw ExceptionExternal(problem);
#endif

#if defined(PETSC_Z)
    if(m_factorised) {
      // already factorised -- so delete and re-create below
      cleanup();
    }

    // store a boolean to indicate that we
    m_factorised = true;
    PetscInt Istart,Iend,n;
    Mat A;

    // size of the (assumed square) matrix
    n = m_pA -> nrows();
    /*
       Create parallel vectors.
    */
    VecCreate(PETSC_COMM_WORLD,&m_petsc_B);
    VecSetSizes(m_petsc_B,PETSC_DECIDE,m_pA->nrows());
    VecSetFromOptions(m_petsc_B);
    VecDuplicate(m_petsc_B,&m_petsc_x);

    // configure the A matrix
    MatCreate(PETSC_COMM_WORLD,&A);
    // set A to be an nxn matrix
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetFromOptions(A);

    // get: all_rows_nnz[i] is the number of nonzero elts in row i
    PetscInt* all_rows_nnz = new PetscInt[ n ];
    m_pA -> nelts_all_rows(all_rows_nnz);

    // pre-allocate memory using the number of non-zero elts
    // in each row (the 0 is ignored here)
    MatSeqAIJSetPreallocation(A, 0, all_rows_nnz);
    MatSetUp(A);

    /*
       Currently, all PETSc parallel matrix formats are partitioned by
       contiguous chunks of rows across the processors.  Determine which
       rows of the matrix are locally owned.
    */
    MatGetOwnershipRange(A,&Istart,&Iend);
    // populate the A matrix using the CppNoddy sparse matrix data
    for(PetscInt i = Istart; i<Iend; ++i) {
      // move the matrix data into PETSc format 1 row at a time
      std::size_t nelts_in_row = all_rows_nnz[i];
      // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
      PetscInt* cols = new PetscInt[nelts_in_row];
      // store the non-zero elts in this row
      PetscScalar* storage = new PetscScalar[nelts_in_row];
      // get the data from the CppNoddy sparse matrix structure
      m_pA -> get_row_petsc(i, storage, cols);
      MatSetValues(A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
      // delete temp storage made in the conversion
      delete[] cols;
      delete[] storage;
    }

    /*
       Assemble matrix, using the 2-step process:
         MatAssemblyBegin(), MatAssemblyEnd()
       Computations can be done while messages are in transition
       by placing code between these two statements.
    */
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  Create the linear solver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    */
    KSPCreate(PETSC_COMM_WORLD,&m_petsc_ksp);
    KSPSetOperators(m_petsc_ksp,A,A);
    KSPSetType(m_petsc_ksp,KSPPREONLY);
    PetscInt  ival,icntl;
    PetscReal val;
    KSPGetPC(m_petsc_ksp,&m_petsc_pc);
    // hardwire a DIRECT SOLVER via MUMPS
    PCSetType(m_petsc_pc,PCLU);
    PCFactorSetMatSolverType(m_petsc_pc,MATSOLVERMUMPS);
    PCFactorSetUpMatSolverType(m_petsc_pc);
    //PCFactorSetMatSolverPackage(m_petsc_pc,MATSOLVERMUMPS);
    //PCFactorSetUpMatSolverPackage(m_petsc_pc);
    /* call MatGetFactor() to create F */
    PCFactorGetMatrix(m_petsc_pc,&m_petsc_F);

    /* sequential ordering */
    //icntl = 7;
    //ival = 2;
    //MatMumpsSetIcntl(m_petsc_F,icntl,ival);

    /* threshhold for row pivot detection */
    //MatMumpsSetIcntl(m_petsc_F,24,1);
    //icntl = 3;
    //val = 1.e-6;
    //MatMumpsSetCntl(m_petsc_F,icntl,val);

    /* compute determinant of A */
    // MatMumpsSetIcntl(m_petsc_F,33,1);
    /* not used unless we initialise PETSc using the command line options */
    // KSPSetFromOptions(m_petsc_ksp);

    /* Get info from matrix factors */
    KSPSetUp(m_petsc_ksp);

    MatDestroy(&A);
    delete[] all_rows_nnz;

#endif // check for PETSC_D/Z
  }




  template <>
  void SparseLinearSystem<std::complex<double> >::solve_using_factorisation() {
#if !defined(PETSC_Z)
    std::string problem;
    problem = "CppNoddy is linked against the DOUBLE version of PETSc\n";
    problem += "but you are trying to solve e a COMPLEX matrix.\n";
    problem += "Recompile with $PETSC_ARCH\n";
    problem += "pointing to a COMPLEX version of the PETSc code.";
    throw ExceptionExternal(problem);
#endif
#if defined(PETSC_Z)
    // size of the (assumed square) matrix
    PetscInt n = m_pA -> nrows();

    // populate the RHS vector using the CppNoddy DenseVector content
    for(PetscInt i = 0; i < n; ++i) {
      VecSetValue(m_petsc_B,i,m_pB->operator[](i),INSERT_VALUES);
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the linear system
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    KSPSolve(m_petsc_ksp,m_petsc_B,m_petsc_x);

    /* We can now gather the parallel result back to ALL processes
      This is temporary as the SparseMatrix is stored on each processes
      and is too dumb for "proper" parallelization */
    Vec y;
    // a scatter context
    VecScatter ctx = 0;
    // map all elts of the parallel vector to a sequential copy
    VecScatterCreateToAll(m_petsc_x,&ctx,&y);
    // scatter it
    VecScatterBegin(ctx,m_petsc_x,y,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_petsc_x,y,INSERT_VALUES,SCATTER_FORWARD);
    // clean up
    VecScatterDestroy(&ctx);
    // this array is a pointer not a copy
    PetscScalar* array;
    VecGetArray(y,&array);
    // now copy to the CppNoddy densevctor
    for(PetscInt i=0; i<n; i++) {
      m_pB -> operator[](i) = array[i];
    }
    // follow the docs and Restore after get
    VecRestoreArray(m_petsc_x,&array);
    VecDestroy(&y);
#endif
  }

  template class SparseLinearSystem<D_complex>
  ;
  template class SparseLinearSystem<double>
  ;

} // end namespace
