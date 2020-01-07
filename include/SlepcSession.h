
#if defined(SLEPC)

#ifndef SLEPCSESSION_H
#define SLEPCSESSION_H

#include <slepc.h>
#include <petsc.h>

namespace CppNoddy {
  
  /* This initialises SLEPc and then finalizes everything on destruction.
     Instantiate this object first, then it will only finalize on
     exit of main() */
  class SlepcSession {

  private:
    SlepcSession(){
    };

  public:

    static SlepcSession* getInstance(int argc, char *argv[]) {
      static SlepcSession instance;
      SlepcInitialize(&argc,&argv,(char*)0,(char*)0);
#if defined(DEBUG)
      PetscPrintf(PETSC_COMM_WORLD, "[DEBUG] Starting a SLEPc session using command line options.\n");
#endif
      return &instance;
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

