/// \file DistributedVector.h
/// A class that constructs a SPARSE/DISTRIBUTED vector using PETSc

#if defined(PETSC_D) || defined(PETSC_Z)

#ifndef DISTRIBUTEDVECTOR_H
#define DISTRIBUTEDVECTOR_H

#include <Exceptions.h>
#include <DenseVector.h>
#include <PetscSession.h>
#include "petsc.h"


namespace CppNoddy {

  /// A class that constructs a SPARSE/DISTRIBUTED vector 
  template <typename _Type>
  class DistributedVector {

   public:

    /// Construct with a set number of elements
    /// \param length The number of elements in the matrix
    DistributedVector(const PetscInt& length) {
#if defined(PARANOID)
      int flag(0);
      MPI_Initialized(&flag);
      if(flag != 1) {
        std::string problem;
        problem = "DistributedVector<> needs PETSc and therefore you must call \n";
        problem += "PetscSession before instantiating the object.\n";
        throw ExceptionRuntime(problem);
      }
#endif
      // create and size a vector m_B
      VecCreate(PETSC_COMM_WORLD,&m_B);
      VecSetSizes(m_B,PETSC_DECIDE,length);
      // add any command line configuration
      VecSetFromOptions(m_B);

      // store the rank/size in the object
      MPI_Comm_size(PETSC_COMM_WORLD,&m_size);
      MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

#if defined(DEBUG)
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Creating a distributed vector\n");
      PetscInt start, end;
      VecGetOwnershipRange(m_B,&start,&end);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Rank = %D Start = %D End = %D \n", m_rank, start, end );
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
#endif
    }

    /// D-tor, make sure we clean up PETSc objects on exit
    ~DistributedVector() {
      // have to destroy on all processes (I think?!)
      VecDestroy(&m_B);
    }

    /// Copy constructor
    DistributedVector(DistributedVector<_Type>& source){
      // VecDuplicate just copies distribution and memory allocation
      VecDuplicate(source.m_B, &m_B);
      m_rank = source.m_rank;
      m_size = source.m_size;
      //
      // VecCopy copies data
      VecCopy(source.m_B, m_B);      
    }

    std::size_t size() const {
      PetscInt nnz(0);
      VecGetSize(m_B, &nnz);
      return std::size_t(nnz);
    }
    
    /// \return 1-norm of the distributed vector
    PetscReal one_norm() const {
      PetscReal norm(0.0);
      VecNorm(m_B,NORM_1,&norm);
      return norm;
    }

    /// \return 2-norm of the distributed vector
    PetscReal two_norm() const {
      PetscReal norm(0.0);
      VecNorm(m_B,NORM_2,&norm);
      return norm;
    }

    /// \return infinity-norm of the distributed vector
    PetscReal inf_norm() const {
      PetscReal norm(0.0);
      VecNorm(m_B,NORM_INFINITY,&norm);
      return norm;
    }

    /// \return A pointer to the PETSc Vector
    Vec* get_pVec(){
      return &m_B;
    }     

    /// \param i The index of the element to set
    /// \param value The value to be put into element i
    void operator()(const PetscInt& i, const PetscScalar& value) {
      PetscInt start, end;
      VecGetOwnershipRange(m_B,&start,&end);
      if ( ( i >= start ) && ( i < end ) ) {
        VecSetValues(m_B,1,&i,&value,INSERT_VALUES);
      } else {
      }
    }

    /// \param i The index of the element to set
    /// \param value The value to be put into element i
    void set_elt(const PetscInt& i, const PetscScalar& value) {
      PetscInt start, end;
      VecGetOwnershipRange(m_B,&start,&end);
      if ( ( i >= start ) && ( i < end ) ) {
        VecSetValues(m_B,1,&i,&value,INSERT_VALUES); 
      } else {
      }
    }
    
    /// \param i The vector of elements to set
    /// \param value The vector of values to be put in the above indices
    void set(const DenseVector<int>& elts, const DenseVector<_Type>& values) {
      PetscInt start, end;
      VecGetOwnershipRange(m_B,&start,&end);
      PetscInt nnz_elts = elts.size();
      VecSetValues(m_B,nnz_elts,&elts[0],&values[0],INSERT_VALUES);
    }

    /// Assemble the matrix
    void final_assembly() {
      VecAssemblyBegin(m_B);
      VecAssemblyEnd(m_B);
    }

    /// View the matrix on stdout
    void view() {
      PetscInt start, end;
      VecGetOwnershipRange(m_B,&start,&end);
      PetscPrintf(PETSC_COMM_WORLD, "Vector view to stdout:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "Processor rank = %D, start = %D end = %D \n", m_rank, start, end-1 );
      PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
      VecView(m_B,PETSC_VIEWER_STDOUT_WORLD);
    }
    
  private:

    // PETSc matrix
    Vec m_B;

    // processor rank
    PetscMPIInt m_rank, m_size;
    
  }
  ; // END CLASS
} // end namespace

#endif // include guard

#endif // check for PETSC_D or PETSC_Z
