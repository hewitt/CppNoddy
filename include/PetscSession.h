
#if defined(PETSC_D) || defined(PETSC_Z)

#ifndef PETSCSESSION_H
#define PETSCSESSION_H

#include <string.h>
#include "petsc.h"

namespace CppNoddy {


  /* This initialises PETSc and then finalizes everything on destruction.
     Instantiate this object first, then it will only finalize on
     exit of main() */
  class PetscSession {

  private:
    /// Constructor is private -- there is only 1 instance
    PetscSession() {
    };
    
  public:

    static PetscSession* getInstance(int argc, char *argv[]) {
      static PetscSession instance;   
      PetscInitialize(&argc,&argv,(char*)0,(char*)0);
#if defined(DEBUG)
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Starting a PETSc session using command line options.\n");
#endif
      return &instance;
    }

     ~PetscSession() {
#if defined(DEBUG)
       PetscSynchronizedPrintf(PETSC_COMM_WORLD, "[DEBUG] Exiting the PETSc session.\n");
#endif
       PetscFinalize();
     }

    
    
  };

} // end CppNoddy namespace

#endif // include guard

#endif // check for PETSC_D and PETSC_Z

