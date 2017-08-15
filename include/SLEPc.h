/// This is a singleton instance used to initialize the SLEPc library.
/// Repeated SLEPCinitialize is OK, but not SLEPcFinalize.
#ifndef SLEPC_H
#define SLEPC_H

#include <slepceps.h>

class SLEPc
{
private:
  static bool instanceFlag;
  static SLEPc *single;
  MPI_Comm COMM;
  SLEPc()
  {
    // Singleton : private constructor
    // initialize the SLEPc library
    SlepcInitialize(0,NULL,NULL,NULL);
    COMM = MPI_COMM_WORLD;
  }
    
public:
    static SLEPc* getInstance();
    
    MPI_Comm get_Comm() const;
    
    ~SLEPc()
    {
        instanceFlag = false;
        SlepcFinalize();
    }
};

bool SLEPc::instanceFlag = false;

SLEPc* SLEPc::single = NULL;

SLEPc* SLEPc::getInstance()
{
    if( !instanceFlag )
    {
        single = new SLEPc();
        instanceFlag = true;
        return single;
    }
    else
    {
        return single;
    }
}

MPI_Comm SLEPc::get_Comm() const
{
  return COMM;
}

#endif


