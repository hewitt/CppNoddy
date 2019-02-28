/// \file ArcLength_base.h
/// A base class for arclength-capable solvers. This defines
/// all the usual get/set methods for arclength parameters and the
/// additional residual used when augmenting the system to allow
/// for iteration on the arclength.

#ifndef ArcLength_base_H
#define ArcLength_base_H

#include <Uncopyable.h>
#include <DenseVector.h>

namespace CppNoddy {

  template <typename _Type>
  class ArcLength_base : private Uncopyable {

   public:

    ArcLength_base() :
      ARCSTEP_MULTIPLIER(2.0),
      INITIALISED(false),
      THETA(0.5),
      DESIRED_ARC_PROPORTION(0.5),
      RESCALE_THETA(false)
    {}

    virtual ~ArcLength_base()
    {}

    /// Initialise the class ready for arc-length continuation.
    /// This will ensure that x is a solution at the current
    /// parameter value and compute the required derivatives
    /// with respect to the parameter.
    /// \param x The initial guess at a solution
    /// \param p A pointer to the required parameter
    /// \param length Initial arc-length step to take
    /// \param max_length Maximum absolute arc length step permitted
    void init_arc(DenseVector<_Type> x,
                  _Type* p,
                  const double& length,
                  const double& max_length);

    /// Compute a solution for that state & parameter variables
    /// that are an arc-length 'ds' from the current state.
    ///virtual void arclength_solve( DenseVector<_Type> &x ) = 0;

    /// Solve the system for a given initial guess.
    /// \param x The initial guess for the system.
    virtual void solve(DenseVector<_Type> &x) = 0;

    /// Return a handle to the arclength step. This is normally done via the init_arc
    /// method, but one can use this  to change the setting
    /// post-initialisation.
    double& ds();

    /// Used to set the multiplication constant used when increasing
    /// or decreasing the arclength step.
    double& arcstep_multiplier();

    /// Handle to the RESCALE_THETA flag. If true the arclength theta
    /// will be rescaled to balance the relative importance of
    /// the state variables and the parameter.
    bool& rescale_theta();

    /// Set the arclength theta parameter.
    double& theta();

    /// Handle to the desired proportion of the parameter to be
    /// used in the arc length solver. If this parameter is 1,
    /// then the system is effectively doing zero-order continuation
    /// in the secondary parameter.
    double& desired_arc_proportion();

   protected:

    /// A method called by arclength_solve and init_arc
    /// which stores the current converged state and parameter
    /// and hence computes the derivatives w.r.t the arc-length.
    void update(const DenseVector<_Type>& x);

    /// The extra constraint that is to be used to replace the
    /// unknown arc length
    /// \return The residual of this constraint
    double arclength_residual(const DenseVector<_Type>& x) const;

    /// The derivative of the arclength_residual function with respect
    /// to each of the state variables. When arc-length continuing boundary
    /// value problems you would NOT want to finite-difference this
    /// since vector lengths will be of the order of 10^3.
    /// \param x The state vector.
    /// \return The derivative of the arclength_residual function with
    /// respect to each state variable.
    DenseVector<_Type> Jac_arclength_residual(DenseVector<_Type>& x) const {
      // a by-hand differentiation of 'arclength_residual' provides
      // the follwoing result.
      DenseVector<_Type> Jx(x - LAST_X);
      Jx *= THETA / (x.size() * (x - LAST_X).two_norm());
      return Jx;
    }

    /// Automatically update the Keller THETA such that the
    /// proportion of the arclength obtained from the parameter
    /// is the desired value. This method will only have any
    /// effect of the RESCALE_THETA flag is true.
    void update_theta(const DenseVector<_Type>& x);

    /// pointer to the parameter in arclength solves
    _Type *p_PARAM;
    /// state variable at the last computed solution
    DenseVector<_Type> LAST_X;
    /// derivative of the state variable w.r.t. arc length
    DenseVector<_Type> X_DERIV_S;
    /// parameter value at the last computed solution
    _Type LAST_PARAM;
    /// derivative of the parameter w.r.t arc length
    _Type PARAM_DERIV_S;
    /// size of the arc length step
    double DS;
    /// maximum arc length step to be taken
    double MAX_DS;
    /// step change multiplier
    double ARCSTEP_MULTIPLIER;
    /// for the arc-length solver - to show it has been initialised
    bool INITIALISED;
    /// the arclength theta
    double THETA;

   private:

    /// the desired proportion of the arc length for the parameter
    double DESIRED_ARC_PROPORTION;
    /// If set, then the arclength theta will be rescaled at each step
    bool RESCALE_THETA;
  };

}

#endif // ArcLength_base_H

