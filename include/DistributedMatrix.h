/// \file DistributedMatrix.h
/// A matrix class that constructs a SPARSE/DISTRIBUTED matrix 
/// using PETSc

#if defined(PETSC_D) || defined(PETSC_Z)

#ifndef DISTRIBUTEDMATRIX_H
#define DISTRIBUTEDMATRIX_H

#include <Exceptions.h>
#include <DenseVector.h>
#include <SparseVector.h>
#include "petsc.h"


namespace CppNoddy {

  /// A matrix class that constructs a SPARSE matrix as a
  /// row major std::vector of SparseVectors.
  template <typename _Type>
  class DistributedMatrix {

   public:

    /// Construct with a set number of rows
    /// \param rows The number of rows in the matrix
    /// \param cols The number of columns in the matrix
    /// \param nd   The number of diagonal entries per row
    /// \param od   The number of off-diagonal entries per row
    DistributedMatrix(const PetscInt& rows, const PetscInt& cols,
                     const PetscInt& nd, const PetscInt& od ) {
      /* "diagonal" vs "off-diagonal" is as defined by PETSc documentation
         -- it depends on the division of the matrix into processor units
         Any columns outside row_start and row_end are OFF diagonal.
      */
#if defined(PETSC_D) || defined(PETSC_Z)
      int flag(0);
      MPI_Initialized(&flag);
      if(flag != 1) {
        std::string problem;
        problem = "DistributedMatrix<> needs PETSc and therefore you must call \n";
        problem += "PetscInitialize before instantiating the object.\n";
        throw ExceptionRuntime(problem);
      }
#else  
      std::string problem;
      problem = "DistributedMatrix<> needs PETSc.\n";
      throw ExceptionRuntime(problem);
#endif

      // temp storage for a single row
      m_temp_row_storage = SparseVector<_Type> (nd+od);

      // set A to be an rows x cols matrix
      MatCreate(PETSC_COMM_WORLD,&m_A);
      MatSetType(m_A,MATMPIAIJ);
      MatSetSizes(m_A,PETSC_DECIDE,PETSC_DECIDE,rows,cols);
      /* getting the preallocation "right" is key to decent performance */
      MatMPIAIJSetPreallocation(m_A, nd, NULL, od, NULL);
      // add any command line configuration
      MatSetFromOptions(m_A);
      MatSetUp(m_A);

      // store the rank/size in the object
      MPI_Comm_size(MPI_COMM_WORLD,&m_size);
      MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);

#if defined(DEBUG)
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Creating a distributed matrix\n");
      PetscInt row_start, row_end;
      MatGetOwnershipRange(m_A,&row_start,&row_end);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Rank = %D Start_row = %D End_row = %D \n", m_rank, row_start, row_end-1 );
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
#endif
    }

    /// D-tor, make sure we clean up PETSc objects on exit
    ~DistributedMatrix(){
      // have to destroy on all processes (I think?!)
      MatDestroy(&m_A);
    }

    /// \return A pointer to the PETSc Mat member data
    Mat* get_pMat(){
      return &m_A;
    }

    void blank() {
      MatZeroEntries(m_A);
    }

    /// \param row The row of the element to set
    /// \param col The column of the element to set
    /// \param value The value to be put into element (row,column)
    void operator()(const PetscInt& row, const PetscInt& col, const PetscScalar& value) {      
      set_elt( row, col, value );
    }

    _Type& operator()(const PetscInt& row, const PetscInt& col) {
      return m_temp_row_storage[i];
    }
    
    /// \param row The row of the element to set
    /// \param col The column of the element to set
    /// \param value The value to be put into element (row,column)
    void set_elt(const PetscInt& row, const PetscInt& col, const PetscScalar& value) {
      PetscInt row_start, row_end;
      MatGetOwnershipRange(m_A,&row_start,&row_end);
      if ( ( row >= row_start ) && ( row < row_end ) ) {
        //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Processor rank = %D, row = %D col = %D \n", m_rank, row, col );
        MatSetValues(m_A,1,&row,1,&col,&value,INSERT_VALUES);
      } else {
      }
    }
    
    /// \param row The row of the element to set
    /// \param cols The vector of column indices to set
    /// \param value The vector of values to be put in the above columns
    void set_row(const PetscInt& row, const DenseVector<int>& cols, const DenseVector<_Type>& values) {
      PetscInt row_start, row_end;
      MatGetOwnershipRange(m_A,&row_start,&row_end);
      if ( ( row >= row_start ) && ( row < row_end ) ) {
        PetscInt nnz_cols = cols.size();
        //PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Processor rank = %D, row = %D cols = %D to %D \n", m_rank, row, cols[0], nnz_cols );
        MatSetValues(m_A,1,&row,nnz_cols,&cols[0],&values[0],INSERT_VALUES);
      } else {
      }
      //PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    }

    /// Assemble the matrix
    void final_assembly() {
      MatAssemblyBegin(m_A,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(m_A,MAT_FINAL_ASSEMBLY);
    }

    /// View the matrix on stdout
    void view() {
      PetscInt row_start, row_end;
      MatGetOwnershipRange(m_A,&row_start,&row_end);
      PetscPrintf(PETSC_COMM_WORLD, "Matrix view to stdout:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Processor rank = %D, start row = %D end row = %D \n", m_rank, row_start, row_end-1 );
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, " Diagonal entry is a %Dx%D square.\n",row_end-row_start,row_end-row_start);
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
      MatView(m_A,PETSC_VIEWER_STDOUT_WORLD);
    }

    
  private:
    // temp storage for building a row at a time
    SparseVector<_Type> m_temp_row_storage;
    
    // PETSc matrix
    Mat m_A;

    // processor rank
    PetscMPIInt m_rank, m_size;
  }
  ; // END CLASS
} // end namespace

#endif // include guard

#endif // check for PETSC_D or PETSC_Z
