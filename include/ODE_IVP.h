/// \file ODE_IVP.h
/// A templated class for an \f$n\f$-th order system
/// of ODEs that form the IVP:
/// \f[ \underline{ \dot f} (t) = \underline R ( \underline f (t), t) \,, \f]
/// where \f$ \underline f (0) \f$ is known.

#ifndef ODE_IVP_H
#define ODE_IVP_H

#include <Types.h>
#include <Equation.h>
#include <OneD_Node_Mesh.h>
#include <Uncopyable.h>

namespace CppNoddy
{

  /// A templated object for real/complex vector system
  /// of first-order ordinary differential equations.

  template <typename _Type>
  class ODE_IVP : private Uncopyable
  {
  public:

    /// The class is defined by a vector function for the system.
    /// \param equation_ptr A pointer to an inherited Equation object.
    /// \param x_init The starting point of the domain for the ODE.
    /// \param x_final The end point of the domain for the ODE.
    /// \param num_of_steps A maximum/default number of steps to be taken.
    ODE_IVP( Equation<_Type > *equation_ptr,
             const double &x_init,
             const double &x_final,
             const std::size_t &num_of_steps );

    ~ODE_IVP();

    /// A fixed step 4th order Runge-Kutta method.
    /// \param u An NVector of initial values.
    /// \return An NVector of the final values.
    DenseVector<_Type> shoot4( DenseVector<_Type> u );

    /// A Runge-Kutta-Fehlberg integrator.
    /// \param u An Nvector of initial values.
    /// \param tol The tolerance used in choosing the step length.
    /// \param h_init The initial step length.
    /// \return An NVector of the final values.
    DenseVector<_Type> shoot45( DenseVector<_Type> u, const double& tol, const double& h_init );

    /// Get the number of steps taken.
    /// \return The number of steps in the integration
    unsigned get_count() const;

    /// Return the history of the stepped solution.
    /// \return A reference to the mesh solution
    OneD_Node_Mesh<_Type>& get_mesh();

    /// Set the flag that determines if the solution 'history'
    /// should be stored.
    /// \param flag The boolean flag to set
    void set_store_soln( bool flag );

  private:

    /// Number of steps taken in the integrator methods
    unsigned COUNT;
    /// initial and final step points
    double X_INIT, X_FINAL;
    /// initial step size
    double H_INIT;
    /// maximum number of steps to take
    std::size_t N;
    /// The function associated with this instance.
    Equation<_Type > *p_EQUATION;
    /// In-class storage of the solution
    OneD_Node_Mesh<_Type> SOLN;
    /// Should we store the whole solution?
    bool STORE_SOLN;
  }
  ; // end class


} // end namespace


#endif

