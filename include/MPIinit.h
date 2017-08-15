/// This is a singleton instance used to initialize the MPI library.
/// Repeated MPIinitialize is OK, but not MPIFinalize.
#ifndef MPIinit_H
#define MPIinit_H

#ifdef INC_MPI

#include <mpi.h>

class MPIinit
{
private:
  static bool instanceFlag;
  static MPIinit *single;
  MPI_Comm COMM;
  MPIinit()
  {
    // Singleton : private constructor
    // initialize the SLEPc library
    MPI_Init(NULL, NULL);
    COMM = MPI_COMM_WORLD;
  }
    
public:
    static MPIinit* getInstance();
    
    MPI_Comm get_Comm() const;
    
    ~MPIinit()
    {
        instanceFlag = false;
        MPI_Finalize();
    }
    
};

bool MPIinit::instanceFlag = false;

MPIinit* MPIinit::single = NULL;

MPIinit* MPIinit::getInstance()
{
    if( !instanceFlag )
    {
        single = new MPIinit();
        instanceFlag = true;
        return single;
    }
    else
    {
        return single;
    }
}

MPI_Comm MPIinit::get_Comm() const
{
  return COMM;
}

#endif // INC_MPI


#endif // MPIinit_H
