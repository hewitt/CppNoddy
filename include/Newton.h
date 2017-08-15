/// \file Newton.h
/// A vector NEWTON iteration class.
/// This allows for Newton iteration to be performed
/// for a vector function of a vector unknown.
/// Use templates to allow double or complex.

#ifndef NEWTON_H
#define NEWTON_H

#include <complex>

#include <DenseMatrix.h>
#include <DenseVector.h>
#include <Residual.h>
#include <ArcLength_base.h>
#include <DenseLinearSystem.h>
#include <Utility.h>

namespace CppNoddy
{
  /// A vector NEWTON iteration class.
  /// This allows for Newton iteration to be performed
  /// for a vector function of a vector unknown.
  /// Use templates to allow double or complex.
  template <typename _Type>
  class Newton : public ArcLength_base< _Type >
  {

  public:

    /// Constructor
    /// \param ptr_to_residual_object A pointer to an inherited Residual object
    /// \param max_steps The maximum number of iteration steps.
    /// \param tolerance A tolerance used as a convergence criterion.
    /// \param derivative_step A step length used to compute derivatives.
    explicit Newton( Residual<_Type >* ptr_to_residual_object,
                     unsigned max_steps = 20,
                     double tolerance = 1.e-8,
                     double derivative_step = 1.e-8 );

    /// The Newton iteration method.
    /// \param x An initial guess vector and returns the solution via this too.
    void iterate( DenseVector<_Type>& x );

    /// If set then the system will monitor the sign of determinant of the
    /// Jacobian matrix and cause an ExceptionBifurcation when it changes
    /// sign.
    /// \param flag The value to be set.
    void set_monitor_det( bool flag )
    {
      MONITOR_DET = flag;
    }

    /// Solve the system for an initial guess by Newton iteration, this
    /// method is inherited from the ArcLength_base class and this simply
    /// points it to the iteration method.
    /// \param x An initial guess vector, the solution overwrites it.
    void solve( DenseVector<_Type>& x )
    {
      iterate( x );
    }

    /// Arc-length solve the system. Before this can be called the
    /// arc_init method should have been called in order to ensure
    /// we know a solution and have derivatives w.r.t. the arc-length
    /// parameter.
    void arclength_solve( DenseVector<_Type>& x );

  private:
    /// tolerance in the iteration convergence test
    double TOL;
    /// a derivative step
    double DELTA;
    /// maximum number of iterations to be taken
    unsigned MAX_STEPS;
    /// last sign of determinant of the Jacobian
    int LAST_DET_SIGN;
    /// pointer to the residual object
    Residual<_Type >* p_RESIDUAL;
    /// flag set to monitor determinant of the Jacobian
    bool MONITOR_DET;
  };


}
#endif // NEWTONV_H
