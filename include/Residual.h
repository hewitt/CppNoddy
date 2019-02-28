/// \file Residual.h
/// A specification of a (double/complex) VECTOR residual class. To be
/// inherited by anything that is to be passed to the Newton object.

#ifndef RESIDUAL_H
#define RESIDUAL_H

#include <Timer.h>
#include <DenseVector.h>
#include <DenseMatrix.h>

namespace CppNoddy {
  /// A base class to be inherited by objects that define residuals
  template <typename _Type>
  class Residual {
   public:
    /// Constructor for a 'square' residual object
    /// that is, N residuals for N unknowns.
    /// \param order The order of the residual vector
    Residual(const unsigned& order);

    /// Constructor for a 'non-square' residual object
    /// that is, there are less residual constraints than unknowns.
    /// \param order The number of residuals
    /// \param nvars The number of unknowns/variables
    Residual(const unsigned& order, const unsigned& nvars);

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Residual();

    /// Update the Residual object for the current set of state variables
    /// \param state The state at which to set the residual object
    void update(const DenseVector<_Type>& state);

    /// Return a handle to the residuals corresponding
    /// to the last update state
    const DenseVector<_Type>& residual() const;

    /// Retrun a handle to the Jacobian of the residual
    /// corresponding to the last update state
    const DenseMatrix<_Type>& jacobian() const;

    /// \return A handle to the step size used when finite-differencing
    /// the residual vector to obtain the Jacobian
    _Type& delta();

    /// \return A handle to the step size used when finite-differencing
    /// the residual vector to obtain the Jacobian
    const _Type& delta() const;

    /// Get the order of the residual vector
    /// \return The order/length of the residual vector
    unsigned get_order() const;

    /// Get the number of variables that this residual condition
    /// is defined for.
    /// \return The number of variables for the residual
    unsigned get_number_of_vars() const;

    /// A blank virtual residual function method.
    /// \param state The unknown variable.
    /// \param f The residual function f(x).
    virtual void residual_fn(const DenseVector<_Type>& state, DenseVector<_Type>& f) const {
      std::string problem;
      problem = "The Residual::residual_fn method has not been implemented.\n";
      problem += "You have to implement this method to define the residual.\n";
      throw ExceptionRuntime(problem);
    }

   protected:
    /// Because the residual evaluation at the current state is assumed
    /// to have already been done by the 'update' method, this routine is
    /// protected. This default uses a finite-differenced Jacobian.
    /// You can overload this to provide an analytic Jacobian if you wish
    /// \param state Dummy state vector ... this is only here to make overloading
    ///    easier, the 'last_x' member data is assumed in the default FD method.
    /// \param jac The NxN matrix Jacobian where the equation_fn
    ///    is a vector function of length N.
    virtual void jacobian(const DenseVector<_Type>& state, DenseMatrix<_Type>& jac) const;

    /// Jacobian for the last state vector
    DenseMatrix<_Type> JAC_AT_LAST_STATE;
    /// Residual for the last state vector
    DenseVector<_Type> FN_AT_LAST_STATE;
    /// The last state vector
    DenseVector<_Type> LAST_STATE;
    /// A default step for FD computation of the Jacobian
    _Type DELTA;
    /// The order of the system of equations
    unsigned ORDER_OF_SYSTEM;
    /// The number of elements in the state vector
    unsigned NUMBER_OF_VARS;
#ifdef TIME
    Timer T_UPDATER;
#endif
  }
  ;  // end class

  template <typename _Type>
  inline void Residual<_Type>::update(const DenseVector<_Type>& x) {
#ifdef TIME
    T_UPDATER.start();
#endif
    LAST_STATE = x;
    residual_fn(LAST_STATE, FN_AT_LAST_STATE);
    jacobian(LAST_STATE, JAC_AT_LAST_STATE);
#ifdef TIME
    T_UPDATER.stop();
#endif
  }

  template <typename _Type>
  inline const DenseVector<_Type>& Residual<_Type>::residual() const {
    return FN_AT_LAST_STATE;
  }

  template <typename _Type>
  inline const DenseMatrix<_Type>& Residual<_Type>::jacobian() const {
    return JAC_AT_LAST_STATE;
  }

  template <typename _Type>
  inline unsigned Residual<_Type>::get_order() const {
    return ORDER_OF_SYSTEM;
  }

  template <typename _Type>
  inline unsigned Residual<_Type>::get_number_of_vars() const {
    return NUMBER_OF_VARS;
  }

  template <typename _Type>
  inline _Type& Residual<_Type>::delta() {
    return DELTA;
  }

  template <typename _Type>
  inline const _Type& Residual<_Type>::delta() const {
    return DELTA;
  }

}   // end namespace

#endif
