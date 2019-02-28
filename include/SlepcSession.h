
#if defined(SLEPC)

#include <slepc.h>
#include <petsc.h>

#ifndef SLEPCSESSION_H
#define SLEPCSESSION_H

namespace CppNoddy {
  
  /* This initialises SLEPc and then finalizes everything on destruction.
     Instantiate this object first, then it will only finalize on
     exit of main() */
  class SlepcSession {
    /// We could be fancy and make this a singleton object
  public:
    SlepcSession(int argc, char *argv[]){
      SlepcInitialize(&argc,&argv,(char*)0,(char*)0);
#if defined(DEBUG)
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Starting a SLEPc session.\n");
      PetscMPIInt rank, size;
      MPI_Comm_size(PETSC_COMM_WORLD, &size);
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Total number of processors = %d.\n",size);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Processor %d of %d reporting in.\n",rank,size);
#endif
    }

    ~SlepcSession(){
#if defined(DEBUG)
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Exiting the PETSc session.\n");
#endif
      SlepcFinalize();
    }
  };

} // end CppNoddy namespace

#endif // include guard

#endif // check for SLEPC

