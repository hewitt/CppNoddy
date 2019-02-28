/// \file OneD_Hyperbolic_System.h

#ifndef ONED_HYPERBOLIC_SYSTEM_H
#define ONED_HYPERBOLIC_SYSTEM_H

#include <Types.h>
#include <Uncopyable.h>
#include <Timer.h>

namespace CppNoddy {

  /// A class to represent a one dimensional hyperbolic system of equations.
  class OneD_Hyperbolic_System : private Uncopyable {

    typedef std::vector<bool> bool_vec;

   public:

    /// \param order The order of the hyperbolic system
    explicit OneD_Hyperbolic_System(const unsigned& order) : ORDER_OF_SYSTEM(order) {
    }

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~OneD_Hyperbolic_System()
    {}

    /// A virtual flux function
    /// \param x The spatial coordinate
    /// \param q The unknowns
    /// \param f The flux function
    virtual void flux_fn(const double &x, const DenseVector<double> &q, DenseVector<double> &f) const {
      std::string problem;
      problem = "The Hyperbolic_Conservative_System::flux_fn method has not been implemented.\n";
      problem += "You have to implement this method to define the system.\n";
      throw ExceptionRuntime(problem);
    }

    /// A virtual function function to define the Jacobian of the
    /// flux function. The default method uses first-order finite
    /// differencing to compute the Jacobian if not otherwise specified
    /// by the user.
    /// \param x The position
    /// \param q The unknowns
    /// \param J The Jacobian of the flux function
    virtual void Jac_flux_fn(const double &x, const DenseVector<double> &q, DenseMatrix<double> &J) const {
      /// first order differencing is the default unless overloaded
      double delta(1.e-8);
      DenseVector<double> state(q);
      DenseVector<double> temp1(ORDER_OF_SYSTEM, 0.0);
      DenseVector<double> temp2(ORDER_OF_SYSTEM, 0.0);
      flux_fn(x, q, temp1);
      // default is to FD the Jacobian
      for(std::size_t i = 0; i < ORDER_OF_SYSTEM; ++i) {
        state[ i ] += delta;
        flux_fn(x, state, temp2);
        state[ i ] -= delta;
        J.set_col(i, (temp2 - temp1) / delta);
      }
    }

    /// A virtual method that is used to bound the shock speed and
    /// must be implemented by the user.
    /// \param q The unknowns
    /// \return A bound on the maximum wave speed
    virtual double max_charac_speed(const DenseVector<double> &q) const {
      std::string problem;
      problem = "The Hyperbolic_Conservative_System::max_shock_speed method has not\n";
      problem += "been implemented. You have to implement this method to define the system.\n";
      throw ExceptionRuntime(problem);
    }

    /// Define the edge boundary conditions.
    /// \param face_index An index for the face:
    /// -1=left +1=right for OneD_TVDLF_Mesh
    /// \param x The position vector along the face
    /// \param q The unknowns specified along the face
    /// \param t A time for unsteady edge conditions
    virtual bool_vec edge_values(const int face_index, const double& x, DenseVector<double>& q, const double &t = 0.0) const {
      return std::vector<bool>(ORDER_OF_SYSTEM, false);
    }

    virtual void source_fn(const double &x, const DenseVector<double> &q, const DenseVector<double> &slope, DenseVector<double>& r) const {
    }

    unsigned get_order() {
      return ORDER_OF_SYSTEM;
    }

   protected:

    /// The order of the system of equations
    const std::size_t ORDER_OF_SYSTEM;
  }
  ; // end class

} // end namespace

#endif
