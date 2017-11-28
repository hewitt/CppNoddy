/// \file SparseLinearSystem.cpp
/// Implementation for the LinearSystem class

#include <vector>
#include <set>

#include <SparseLinearSystem.h>
#include <Exceptions.h>
#include <Types.h>
#include <Timer.h>

#ifdef INC_MPI
  #include "mpi.h"
#endif

namespace CppNoddy
{

  template <typename _Type>
  SparseLinearSystem<_Type>::SparseLinearSystem( SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which ) : LinearSystem_base(),
      MIN_PIV( 1.e-12 ), factorised_(false)
  {
    p_A = Aptr;
    p_B = Bptr;
    VERSION = which;
    //
    if ( ( VERSION != "petsc" ) && ( VERSION != "native" ) )
    {
      std::cout << "Solver type requested: " << VERSION << "\n";
      std::string problem;
      problem = "The SparseLinearSystem has been instantiated with an unrecognised\n";
      problem += "request for a solver type. Options: 'native','petsc'. \n";
      throw ExceptionRuntime( problem );
    }
    #ifdef INC_MPI
      if ( VERSION == "petsc" )
      {
        int flag(0);
        MPI_Initialized( &flag );
        if ( flag != 1 )
        {
          std::string problem;
          problem = "The SparseLinearSystem has been instantiated for a petsc solver.\n";
          problem += "You must call PetscInitialize before calling the petsc solver.\n";
          throw ExceptionRuntime( problem );
        }
      }
    #endif
    #ifdef INC_MPI
      MPI_Comm_size(MPI_COMM_WORLD,&size_);
      MPI_Comm_size(MPI_COMM_WORLD,&rank_);
      if ( size_ > 1 )
      {
        std::string problem;
        problem = " The SparseLinearSystem object links to PETSc which makes\n";
        problem += " use of MPI, but you probably won't gain anything yet by\n";
        problem += " using >1 processes. The SparseMatrix object is still too\n";
        problem += " dumb, and will be stored in each process.\n";
        throw ExceptionRuntime( problem );
      }
    #endif
  }

  template<typename _Type>
  SparseLinearSystem<_Type>::~SparseLinearSystem()
  {
    cleanup();
  }

  template<typename _Type>
  void SparseLinearSystem<_Type>::cleanup()
  {
    #if defined(PETSC_D) || defined(PETSC_Z)
      // delete objects used in the factorisation?
      KSPDestroy(&ksp_);
      VecDestroy(&x_);
      VecDestroy(&B_);
    #endif
    factorised_ = false;
  }

  template <typename _Type>
  void SparseLinearSystem<_Type>::solve()
  {
    if ( MONITOR_DET )
    {
      std::string problem;
      problem = "You've asked for the LinearSystem object to monitor \n";
      problem += "the determinant for a sparse matrix. This has not been \n";
      problem += "implemented.\n";
      throw ExceptionRuntime( problem );
    }
    if ( "petsc" == VERSION )
    {
      factorise();
      solve_using_factorisation();
    }
    else // we catch incorrect VERSION choices in the ctor
    {
      solve_native();
    }
  }



  template <typename _Type>
  void SparseLinearSystem<_Type>::solve_native()
  {
    const std::size_t Nr = p_A -> nrows();
    // step through rows
    for ( std::size_t l = 0 ; l < Nr - 1 ; ++l )
    {
      // find max index in column 'l' in the range [l, Nr)
      const std::size_t index = p_A -> max_in_col( l , l , Nr );

      // if index of max elt is not the diagonal then swap the rows
      if ( l != index )
      {
        // switch the index row (with maxel) to current position "l"
        p_A -> row_swap( l, index );
        // switch the elts in RHS R-vector
        p_B -> swap( l, index );
      }

      // The diagonal entry is the first :
      // const _Type diag_entry = matrix[ l ].get( l );
      const _Type diag_entry = ( p_A -> MATRIX[ l ] ).begin() -> second;
#ifdef PARANOID
      if ( std::abs( diag_entry ) < MIN_PIV )
      {
        std::string problem( "The pivot in SparseLinearSystem is under the minimum tolerance.\n" );
        throw ExceptionRuntime( problem );
      }
#endif
      // eliminate all entries below
      for ( std::size_t row = l + 1 ; row < Nr ; ++row )
      {
        // subtract rows: R_row = R_row - mult * R_l
        // but optimise out the zero elements to the left
        // cycle through R_l data
        typename std::map< std::size_t, _Type >::const_iterator pos_ro = ( p_A -> MATRIX[ l ] ).begin();
        typename std::map< std::size_t, _Type >::iterator pos_rw = ( p_A -> MATRIX[ row ] ).begin();
        if ( pos_rw -> first <= l )
        {
          // this should be equivalent to 'MATRIX[ row ].get( l ) / diag_entry;"
          const _Type mult = ( pos_rw -> second ) / diag_entry;
          do
          {
            std::size_t index_rw = pos_rw -> first;
            std::size_t index_ro = pos_ro -> first;
            if ( index_rw == index_ro )
            {
              // element in both vectors
              if ( pos_rw -> first == l )
              {
                // this entry should now be eliminated from the matrix
                ( p_A -> MATRIX[ row ] ).erase( pos_rw );
              }
              else
              {
                pos_rw -> second -= ( pos_ro -> second ) * mult;
              }
              ++pos_rw;
              ++pos_ro;
            }
            if ( index_rw > index_ro )
            {
              // element is in X but not 'this'
              ( p_A -> MATRIX[ row ] ).set( index_ro ) = -( pos_ro -> second ) * mult;
              ++pos_ro;
            }
            if ( index_rw < index_ro )
            {
              // element is in 'this' but not X
              ++pos_rw;
            }
          }
          while ( pos_ro != ( p_A -> MATRIX[ l ] ).end() &&
                  pos_rw != ( p_A -> MATRIX[ row ] ).end() );

          if ( pos_ro != ( p_A -> MATRIX[ l ] ).end() )
          {
            // need to finish the X data
            do
            {
              ( p_A -> MATRIX[ row ] ).set( pos_ro -> first ) = -( pos_ro -> second ) * mult;
              ++pos_ro;
            }
            while ( pos_ro != ( p_A -> MATRIX[ l ] ).end() );
          }
          // this is a scalar operation
          p_B -> operator[] ( row ) -= p_B -> operator[] ( l ) * mult;
        }
      }

    } // close l-loop
#ifdef PARANOID
    // check last row for singular matrix
    // The diagonal entry should be the first:
    // const _Type diag_entry = matrix[ Nr - 1 ].get( Nr - 1 );
    const _Type diag_entry = ( p_A -> MATRIX[ Nr - 1 ] ).begin() -> second;
    if ( std::abs( diag_entry ) < MIN_PIV )
    {
      std::string problem( "The pivot in NSMatrix.GJE is under the minimum tolerance.\n" );
      throw ExceptionRuntime( problem );
    }
#endif
    backsub( *p_A, *p_B );
  }



/*
  $PETSC_ARCH points to a DOUBLE implementation routines are below
*/

template<>
void SparseLinearSystem<double>::factorise()
{
#if !defined(PETSC_D)
  std::string problem;
  problem = "CppNoddy is linked against the COMPLEX version of PETSc\n";
  problem += "but you are trying to factorise a DOUBLE matrix. Either\n";
  problem += "redefine your matrix as complex, or recompile with $PETSC_ARCH\n";
  problem += "pointing to a DOUBLE version of the PETSc code.";
  throw ExceptionExternal( problem );
#endif
#if defined(PETSC_D)
  if (factorised_)
  {
    // already factorised -- so delete and re-create below
    cleanup();
  }

  // store a boolean to indicate that we
  factorised_ = true;
  PetscInt Istart,Iend,n;
  Mat A;

  // size of the (assumed square) matrix
  n = p_A -> nrows();
  /*
     Create parallel vectors.
  */
  VecCreate(PETSC_COMM_WORLD,&B_);
  VecSetSizes(B_,PETSC_DECIDE,p_A->nrows());
  VecSetFromOptions(B_);
  VecDuplicate(B_,&x_);

  // configure the A matrix
  MatCreate(PETSC_COMM_WORLD,&A);
  // set A to be an nxn matrix
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
  MatSetFromOptions(A);

  // get: all_rows_nnz[i] is the number of nonzero elts in row i
  PetscInt* all_rows_nnz = new PetscInt[ n ];
  p_A -> nelts_all_rows( all_rows_nnz );

  // pre-allocate memory using the number of non-zero elts in each row (the 0 is ignored here)
  MatSeqAIJSetPreallocation(A, 0, all_rows_nnz );
  // need to allocate for MPI codes too
  // \todo if we every get MPI running, we need to sort out preallocation
  // MatMPIAIJSetPreallocation(A, 800, NULL, 800, NULL);
  //
  // finish the A definition
  MatSetUp(A);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  MatGetOwnershipRange(A,&Istart,&Iend);
  // populate the A matrix using the CppNoddy sparse matrix data
  for ( PetscInt i = Istart; i<Iend; ++i )
  {
    // move the matrix data into PETSc format 1 row at a time
    std::size_t nelts_in_row = all_rows_nnz[i]; //p_A -> nelts_in_row(i);
    // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
    PetscInt* cols = new PetscInt[all_rows_nnz[i]];
    // store the non-zero elts in this row
    PetscScalar* storage = new PetscScalar[all_rows_nnz[i]];
    // get the data from the CppNoddy sparse matrix structure
    p_A -> get_row_petsc( i, storage, cols );
    MatSetValues(A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
    // delete temp storage made in the conversion
    delete[] cols; delete[] storage;
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
  KSPCreate(PETSC_COMM_WORLD,&ksp_);
  KSPSetOperators(ksp_,A,A);
  KSPSetType(ksp_,KSPPREONLY);
  PetscInt  ival,icntl;
  PetscReal val;
  KSPGetPC(ksp_,&pc_);
  // hardwire a DIRECT SOLVER via MUMPS
  PCSetType(pc_,PCLU);
  PCFactorSetMatSolverPackage(pc_,MATSOLVERMUMPS);
  PCFactorSetUpMatSolverPackage(pc_);
  /* call MatGetFactor() to create F */
  PCFactorGetMatrix(pc_,&F_);

  /* sequential ordering */
  icntl = 7; ival = 2;
  MatMumpsSetIcntl(F_,icntl,ival);

  /* threshhold for row pivot detection */
  MatMumpsSetIcntl(F_,24,1);
  icntl = 3; val = 1.e-6;
  MatMumpsSetCntl(F_,icntl,val);

  /* compute determinant of A */
  // MatMumpsSetIcntl(F_,33,1);
  /* not used unless we initialise PETSc using the command line options */
  KSPSetFromOptions(ksp_);

  /* Get info from matrix factors */
  KSPSetUp(ksp_);

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
  MatDestroy(&A);
  delete[] all_rows_nnz;
#endif // check for PETSC_D/Z
}




template <>
void SparseLinearSystem<double>::solve_using_factorisation()
{
#if !defined(PETSC_D)
  std::string problem;
  problem = "CppNoddy is linked against the COMPLEX version of PETSc\n";
  problem += "but you are trying to solve a DOUBLE matrix. Either\n";
  problem += "redefine your matrix as complex, or recompile with $PETSC_ARCH\n";
  problem += "pointing to a DOUBLE version of the PETSc code.";
  throw ExceptionExternal( problem );
#endif
#if defined(PETSC_D)
  // size of the (assumed square) matrix
  PetscInt n = p_A -> nrows();

  // populate the RHS vector using the CppNoddy DenseVector content
  for ( PetscInt i = 0; i < n; ++i )
  {
    VecSetValue(B_,i,p_B->operator[](i),INSERT_VALUES);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  KSPSolve(ksp_,B_,x_);

  /* We can now gather the parallel result back to ALL processes
    This is temporary as the SparseMatrix is stored on each processes
    and is too dumb for "proper" parallelization */
  Vec y;
  // a scatter context
  VecScatter ctx = 0;
  // map all elts of the parallel vector to a sequential copy
  VecScatterCreateToAll(x_,&ctx,&y);
  // scatter it
  VecScatterBegin(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
  // clean up
  VecScatterDestroy(&ctx);
  // this array is a pointer not a copy
  PetscScalar* array;
  VecGetArray(y,&array);
  // now copy to the CppNoddy densevctor
  for (PetscInt i=0; i<n; i++)
  {
    p_B -> operator[](i) = array[i];
  }
  // follow the docs and Restore after get
  VecRestoreArray(x_,&array);
  #endif
}






/*
  $PETSC_ARCH points to a COMPLEX implementation routines are below
*/

template<>
void SparseLinearSystem<std::complex<double> >::factorise()
{
#if !defined(PETSC_Z)
  std::string problem;
  problem = "CppNoddy is linked against the DOUBLE version of PETSc\n";
  problem += "but you are trying to factorise a COMPLEX matrix.\n";
  problem += "Recompile with $PETSC_ARCH\n";
  problem += "pointing to a COMPLEX version of the PETSc code.";
  throw ExceptionExternal( problem );
#endif
#if defined(PETSC_Z)
  if (factorised_)
  {
    // already factorised -- so delete and re-create below
    cleanup();
  }

  // store a boolean to indicate that we
  factorised_ = true;
  PetscInt Istart,Iend,n;
  Mat A;

  // size of the (assumed square) matrix
  n = p_A -> nrows();
  /*
     Create parallel vectors.
  */
  VecCreate(PETSC_COMM_WORLD,&B_);
  VecSetSizes(B_,PETSC_DECIDE,p_A->nrows());
  VecSetFromOptions(B_);
  VecDuplicate(B_,&x_);

  // configure the A matrix
  MatCreate(PETSC_COMM_WORLD,&A);
  // set A to be an nxn matrix
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
  MatSetFromOptions(A);

  // get: all_rows_nnz[i] is the number of nonzero elts in row i
  PetscInt* all_rows_nnz = new PetscInt[ n ];
  p_A -> nelts_all_rows( all_rows_nnz );

  // pre-allocate memory using the number of non-zero elts in each row (the 0 is ignored here)
  MatSeqAIJSetPreallocation(A, 0, all_rows_nnz );
  // need to allocate for MPI codes too
  // \todo if we every get MPI running, we need to sort out preallocation
  // MatMPIAIJSetPreallocation(A, 800, NULL, 800, NULL);
  //
  // finish the A definition
  MatSetUp(A);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  MatGetOwnershipRange(A,&Istart,&Iend);
  // populate the A matrix using the CppNoddy sparse matrix data
  for ( PetscInt i = Istart; i<Iend; ++i )
  {
    // move the matrix data into PETSc format 1 row at a time
    std::size_t nelts_in_row = all_rows_nnz[i];
    // row i has all_rows_nnz[i] elements that are non-zero, so we store their columns
    PetscInt* cols = new PetscInt[nelts_in_row];
    // store the non-zero elts in this row
    PetscScalar* storage = new PetscScalar[nelts_in_row];
    // get the data from the CppNoddy sparse matrix structure
    p_A -> get_row_petsc( i, storage, cols );
    MatSetValues(A,1,&i,nelts_in_row,cols,storage,INSERT_VALUES);
    // delete temp storage made in the conversion
    delete[] cols; delete[] storage;
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
  KSPCreate(PETSC_COMM_WORLD,&ksp_);
  KSPSetOperators(ksp_,A,A);
  KSPSetType(ksp_,KSPPREONLY);
  PetscInt  ival,icntl;
  PetscReal val;
  KSPGetPC(ksp_,&pc_);
  // hardwire a DIRECT SOLVER via MUMPS
  PCSetType(pc_,PCLU);
  PCFactorSetMatSolverPackage(pc_,MATSOLVERMUMPS);
  PCFactorSetUpMatSolverPackage(pc_);
  /* call MatGetFactor() to create F */
  PCFactorGetMatrix(pc_,&F_);

  /* sequential ordering */
  icntl = 7; ival = 2;
  MatMumpsSetIcntl(F_,icntl,ival);

  /* threshhold for row pivot detection */
  MatMumpsSetIcntl(F_,24,1);
  icntl = 3; val = 1.e-6;
  MatMumpsSetCntl(F_,icntl,val);

  /* compute determinant of A */
  // MatMumpsSetIcntl(F_,33,1);
  /* not used unless we initialise PETSc using the command line options */
  KSPSetFromOptions(ksp_);

  /* Get info from matrix factors */
  KSPSetUp(ksp_);

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
  MatDestroy(&A);
  delete[] all_rows_nnz;
#endif // check for PETSC_D/Z
}




template <>
void SparseLinearSystem<std::complex<double> >::solve_using_factorisation()
{
#if !defined(PETSC_Z)
  std::string problem;
  problem = "CppNoddy is linked against the DOUBLE version of PETSc\n";
  problem += "but you are trying to solve e a COMPLEX matrix.\n";
  problem += "Recompile with $PETSC_ARCH\n";
  problem += "pointing to a COMPLEX version of the PETSc code.";
  throw ExceptionExternal( problem );
#endif
#if defined(PETSC_Z)
  // size of the (assumed square) matrix
  PetscInt n = p_A -> nrows();

  // populate the RHS vector using the CppNoddy DenseVector content
  for ( PetscInt i = 0; i < n; ++i )
  {
    VecSetValue(B_,i,p_B->operator[](i),INSERT_VALUES);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  KSPSolve(ksp_,B_,x_);

  /* We can now gather the parallel result back to ALL processes
    This is temporary as the SparseMatrix is stored on each processes
    and is too dumb for "proper" parallelization */
  Vec y;
  // a scatter context
  VecScatter ctx = 0;
  // map all elts of the parallel vector to a sequential copy
  VecScatterCreateToAll(x_,&ctx,&y);
  // scatter it
  VecScatterBegin(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x_,y,INSERT_VALUES,SCATTER_FORWARD);
  // clean up
  VecScatterDestroy(&ctx);
  // this array is a pointer not a copy
  PetscScalar* array;
  VecGetArray(y,&array);
  // now copy to the CppNoddy densevctor
  for (PetscInt i=0; i<n; i++)
  {
    p_B -> operator[](i) = array[i];
  }
  // follow the docs and Restore after get
  VecRestoreArray(x_,&array);
  VecDestroy(&y);
  #endif
}



  template <typename _Type>
  void SparseLinearSystem<_Type >::backsub( SparseMatrix<_Type> &A, DenseVector<_Type> &B ) const
  {
    const std::size_t Nr( B.size() );
    DenseVector<_Type> x( Nr, 0.0 );
    // This line should be equivalent to:
    // x[ Nr - 1 ] = B[ Nr - 1 ] / matrix[ Nr - 1 ].get( Nr - 1 );
    x[ Nr - 1 ] = B[ Nr - 1 ] / ( A.MATRIX[ Nr - 1 ].begin() -> second );
    // Note the unusual (aka hacked up) row < N termination condition.
    // We can't do "row <= 0" with std::size_t.
    for ( std::size_t row = Nr - 2; row < Nr; --row )
    {
      _Type sum = 0.0;
      for ( typename std::map< std::size_t, _Type >::const_iterator pos_ro =
              A.MATRIX[ row ].begin(); pos_ro != A.MATRIX[ row ].end(); ++pos_ro )
      {
        std::size_t index_ro = pos_ro -> first;
        if ( index_ro > row )
        {
          sum += ( pos_ro -> second ) * x[ index_ro ];
        }
      }
      // The diagonal entry is the first since it's triangular
      x[ row ] = ( B[ row ] - sum ) / ( A.MATRIX[ row ].begin() -> second );
    }
    B = x;
  }






  template class SparseLinearSystem<D_complex>
  ;
  template class SparseLinearSystem<double>
  ;

} // end namespace
