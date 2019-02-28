
#if defined(PETSC_D) || defined(PETSC_Z)

#ifndef PETSCSESSION_H
#define PETSCSESSION_H

namespace CppNoddy {
  
  /* This initialises PETSc and then finalizes everything on destruction.
     Instantiate this object first, then it will only finalize on
     exit of main() */
  class PetscSession {
    /// We could be fancy and make this a singleton object
  public:
    PetscSession(int argc, char *argv[]){
      PetscInitialize(&argc,&argv,(char*)0,(char*)0);
#if defined(DEBUG)
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Starting a PETSc session.\n");
      PetscMPIInt rank, size;
      MPI_Comm_size(PETSC_COMM_WORLD, &size);
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Total number of processors = %d.\n",size);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Processor %d of %d reporting in.\n",rank,size);
#endif
    }

    ~PetscSession(){
#if defined(DEBUG)
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Exiting the PETSc session.\n");
#endif
      PetscFinalize();
    }
  };

} // end CppNoddy namespace

#endif // include guard

#endif // check for PETSC_D and PETSC_Z

