/// \file LinearSystem_base.cpp
/// Implementation for the LinearSystem class

#include <set>

#include <LinearSystem_base.h>
#include <Exceptions.h>
#include <Types.h>

namespace CppNoddy
{

  LinearSystem_base::LinearSystem_base() :
      DET_SIGN( 0 ),
      MONITOR_DET( false )
  {}

  LinearSystem_base::~LinearSystem_base()
  {}

  // MISC. get/set

  int LinearSystem_base::get_det_sign() const
  {
    return DET_SIGN;
  }

  void LinearSystem_base::set_monitor_det( bool flag )
  {
    MONITOR_DET = flag;
  }

  void LinearSystem_base::solve()
  {
    std::string problem;
    problem = "The solve method has not been implemented for\n";
    problem += "the container class you are using here.\n";
    throw ExceptionRuntime( problem );
  }


} // end namespace
