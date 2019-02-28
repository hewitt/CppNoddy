/// \file src/Newton.cpp
/// Implementation of the vector NEWTON iteration class.
/// This allows for Newton iteration to be performed
/// for a vector function of a vector unknown.
/// Use templates to allow double or complex.

#include <complex>
#include <string>

#include <Newton.h>
#include <Residual.h>
#include <Exceptions.h>
#include <DenseLinearSystem.h>
#include <ArcLength_base.h>

namespace CppNoddy {

  template <typename _Type>
  Newton<_Type>::Newton(Residual<_Type >* ptr_to_residual_object,
                        unsigned max_steps,
                        double tolerance,
                        double derivative_step)
    : ArcLength_base< _Type >(),
      TOL(tolerance),
      MAX_STEPS(max_steps),
      p_RESIDUAL(ptr_to_residual_object) {
    p_RESIDUAL -> delta() = derivative_step;
    DELTA = derivative_step;
    MONITOR_DET = false;
    LAST_DET_SIGN = 1;
  }


  template <typename _Type>
  void Newton<_Type>::iterate(DenseVector<_Type>& x) {
    // length of the vector
    unsigned N = x.size();
    // Jacobian
    DenseMatrix<_Type> J(N, N, 0.0);
    // NVectors
    DenseVector<_Type> oldFn(N, 0.0), newFn(N, 0.0);
    // linear system object - native because no LAPACK complex at the moment
    DenseLinearSystem<_Type> system(&J, &oldFn, "native");
    system.set_monitor_det(MONITOR_DET);
    // iteration counter
    unsigned itn = 0;
    do {
      // increment the counter
      ++itn;

      // evaluate the residuals and jacobian at the current state
      p_RESIDUAL -> update(x);
      // get the residuals
      oldFn = p_RESIDUAL -> residual();
#ifdef DEBUG
      std::cout << " DEBUG: starting with |Residuals|  = " << oldFn.inf_norm() << "\n";
#endif

      if((std::abs(oldFn.inf_norm()) < TOL) || (itn == MAX_STEPS)) {
        break;
      }
      // retrieve the current jacobian
      J = p_RESIDUAL -> jacobian();

      // linear solver
      // \todo LU interface to LAPACK for complex matrices
      system.solve();

#ifdef DEBUG
      std::cout << " DEBUG: Iteration number    = " << itn << "\n";
      std::cout << " DEBUG: |Newton correction| = " << oldFn.inf_norm() << "\n";
#endif

      // must *subtract* delta
      x -= oldFn;
    } while(true);

    // More the 'MAX_STEPS' iterations currently triggers a failure.
    if(itn == MAX_STEPS) {
      std::string problem;
      problem = " The Newton.iterate method took too many iterations. \n";
      problem += " At the moment, this is set as a failure. \n";
      throw ExceptionItn(problem, itn, oldFn.inf_norm());
    }
    LAST_DET_SIGN = system.get_det_sign();
    if(MONITOR_DET) {
      if(system.get_det_sign() * LAST_DET_SIGN < 0) {
        std::string problem;
        problem = "[ INFO ] : Determinant monitor has changed signs in ODE_BVP.\n";
        problem += "[ INFO ] : Bifurcation detected.\n";
        throw ExceptionBifurcation(problem);
      }
    }
  }


  template <typename _Type>
  void Newton<_Type>::arclength_solve(DenseVector<_Type>& x) {
#ifdef PARANOID
    if(!this -> INITIALISED) {
      std::string problem;
      problem = "The Newton.arclength_solve method has been called, but \n";
      problem += " you haven't called the arc_init method first. This means \n";
      problem += " that a starting solution & derivatives thereof are unknown. \n";
      problem += " Please initialise things appropriately! ";
      throw ExceptionRuntime(problem);
    }
#endif
#ifdef DEBUG
    std::cout.precision(6);
    std::cout << "[ DEBUG ] : Entering arclength_solve of the Newton class with\n";
    std::cout << "[ DEBUG ] : a parameter of " << *(this -> p_PARAM) << "\n";
#endif
    // backup the state/system in case we fail
    DenseVector<_Type> backup_state(x);
    _Type backup_parameter(*(this -> p_PARAM));

    int det_sign(1);
    // init some local vars
    bool step_succeeded(false);
    unsigned itn(0);
    // make a guess at the next solution
    *(this -> p_PARAM) = this -> LAST_PARAM + this -> PARAM_DERIV_S * this -> DS;
    x = this -> LAST_X + this -> X_DERIV_S * this -> DS;

    // the order of the system
    unsigned N = this -> p_RESIDUAL -> get_order();
    DenseMatrix<_Type> J(N, N, 0.0);
    DenseVector<_Type> R1(N, 0.0);
    DenseVector<_Type> R2(N, 0.0);
    DenseVector<_Type> dR_dp(N, 0.0);
    //
    do {
      // update the residual object to the current guess
      this -> p_RESIDUAL -> update(x);
      // get the Jacobian of the system
      J = this -> p_RESIDUAL -> jacobian();
      R1 = this -> p_RESIDUAL -> residual();
      // get the residual of the EXTRA arclength residual
      double E1 = this -> arclength_residual(x);
      // compute derivatives w.r.t the parameter
      *(this -> p_PARAM) += this -> DELTA;
      this -> p_RESIDUAL -> residual_fn(x, R2);
      double E2 = this -> arclength_residual(x);
      *(this -> p_PARAM)  -= this -> DELTA;
      dR_dp = (R2 - R1) / this -> DELTA;
      _Type dE_dp = (E2 - E1) / this -> DELTA;
      // bordering algorithm
      DenseVector<_Type> y(-R1);
      DenseVector<_Type> z(-dR_dp);
#ifdef LAPACK
      DenseLinearSystem<_Type> sys1(&J, &y, "lapack");
#else
      DenseLinearSystem<_Type> sys1(&J, &y, "native");
#endif
      sys1.set_monitor_det(MONITOR_DET);
      try {
        sys1.solve();
        det_sign = sys1.get_det_sign();
      } catch ( const ExceptionExternal &error ) {
        break;
      }
      // the solve will have overwritten the LHS, so we need to replace it
      J = this -> p_RESIDUAL -> jacobian();
#ifdef LAPACK
      DenseLinearSystem<_Type> sys2(&J, &z, "lapack");
#else
      DenseLinearSystem<_Type> sys2(&J, &z, "native");
#endif
      try {
        sys2.solve();
      } catch( const ExceptionExternal &error ) {
        break;
      }
      DenseVector<_Type> JacE(this -> Jac_arclength_residual(x));
      _Type delta_p = - (E1 + Utility::dot(JacE, y)) /
                      (dE_dp + Utility::dot(JacE, z));
      DenseVector<_Type> delta_x = y + z * delta_p;
      double max_correction = std::max(delta_x.inf_norm(), std::abs(delta_p));
      if(max_correction < this -> TOL) {
        step_succeeded = true;
        break;
      }
      // add the corrections to the state variables
      x += delta_x;
      *(this -> p_PARAM) += delta_p;
      ++itn;
      if(itn > MAX_STEPS) {
        step_succeeded = false;
        break;
      }
    } while(true);

    // is this a successful step?
    if(!step_succeeded) {
#ifdef DEBUG
      std::cout << "[ DEBUG ] : REJECTING STEP \n";
#endif
      // if not a successful step then restore things
      x = backup_state;
      *(this -> p_PARAM) = backup_parameter;
      // restore the residual object to its initial state
      this -> p_RESIDUAL -> update(x);
      // reduce our step length
      this -> DS /= this -> ARCSTEP_MULTIPLIER;
    } else {
      // update the variables needed for arc-length continuation
      this -> update(x);
      if(LAST_DET_SIGN * det_sign < 0) {
        LAST_DET_SIGN = det_sign;
        std::string problem;
        problem = "[ INFO ] : Determinant monitor has changed signs in the Newton class.\n";
        problem += "[ INFO ] : Bifurcation detected.\n";
        throw ExceptionBifurcation(problem);
      } else {
        LAST_DET_SIGN = det_sign;
      }
#ifdef DEBUG
      std::cout << "[ DEBUG ] : Number of iterations = " << itn << "\n";
      std::cout << "[ DEBUG ] : Parameter p = " << *(this -> p_PARAM)
                << "; arclength DS = " << this -> DS << "\n";
#endif
      if(itn >= 7) {
        // converging too slowly, so decrease DS
        this -> DS /= this -> ARCSTEP_MULTIPLIER;
#ifdef DEBUG
        std::cout << "[ DEBUG ] : I decreased DS to " << this -> DS << "\n";
#endif
      }
      if(itn <= 2) {
        if(std::abs(this -> DS * this -> ARCSTEP_MULTIPLIER) < this -> MAX_DS) {
          // converging too quickly, so increase DS
          this -> DS *= this -> ARCSTEP_MULTIPLIER;
#ifdef DEBUG
          std::cout << "[ DEBUG ] : I increased DS to " << this -> DS << "\n";
#endif
        }
      }
    }
  }

  // the templated versions we require are:
  template class Newton<double>
  ;
  template class Newton<std::complex<double> >
  ;

}
