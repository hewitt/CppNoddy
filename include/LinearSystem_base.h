/// \file LinearSystem_base.h
/// Specification of the linear system class.

#ifndef LINEARSYSTEM_BASE_H
#define LINEARSYSTEM_BASE_H

#include <string>

#include <Uncopyable.h>

namespace CppNoddy
{

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  /// This base class stores the get/set methods for
  /// monitoring the determinant of the system and
  /// the tag for selecting the solver.
  class LinearSystem_base : private Uncopyable
  {

  public:

    /// Constructor for a linear system object.
    LinearSystem_base();

    /// Destructor for a linear system object.
    virtual ~LinearSystem_base();

    /// Linear solver
    virtual void solve();

    /// Get the sign of the determinant of the LHS matrix
    /// from the linear system just computed.
    /// \return The sign of the determinant of the
    /// LAST solved system.
    int get_det_sign() const;

    /// Store the sign of the determinant of the LHS matrix
    /// every time a solve is requested on a real system.
    /// \param flag The boolean value to set.
    void set_monitor_det( bool flag );

  protected:

    /// a string ID to pick out the appropriate solver
    std::string VERSION;

    /// the sign of the determinant of the last solved system LHS
    int DET_SIGN;

    /// a flag that determines of the determinant sign should be monitored
    bool MONITOR_DET;

  };

} //end namepsace
#endif
