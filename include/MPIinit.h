/// This is a singleton instance used to initialize the MPI library.
/// Repeated MPIinitialize is OK, but not MPIFinalize.
#ifndef MPIinit_H
#define MPIinit_H

#ifdef INC_MPI

#include <mpi.h>

namespace CppNoddy
{
class MPIinit
{
private:
  static bool instanceFlag;
  static MPIinit *single;
  MPI_Comm COMM;
  MPIinit()
  {
    std::cout << "Initializing MPI\n";
    // Singleton : private constructor
    // initialize the SLEPc library
    MPI_Init(NULL, NULL);
    COMM = MPI_COMM_WORLD;
    //MPI_Comm_dup(MPI_COMM_WORLD,COMM);
  }

public:
    static MPIinit* getInstance();

    MPI_Comm get_Comm() const;

    ~MPIinit()
    {
        std::cout << "Running destructor for MPIinit\n";
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
}

#endif // INC_MPI


#endif // MPIinit_H
