/// \file ODE_BVP.h
/// A specification of a class for an \f$ n^{th} \f$-order ODE BVP defined by
/// \f[ {\underline f}^\prime (y) = {\underline R}( {\underline f}(y), y )\,, \f]
/// subject to \f$ n \f$ Dirichlet conditions defined at \f$ y = y_{left} \f$ or
/// \f$ y_{right} \f$ for some components of \f$ {\underline f}(y) \f$. The system
/// is solved by applying Newton iteration, with the intermediate problem:
/// \f[ {\underline g}^\prime (y) - \frac{\partial {\underline R}}{\partial \underline f} \Big \vert_{\underline F} \,\,{\underline g}(y) = {\underline R}( {\underline F}(y), y) - {\underline F}^\prime (y) \,, \f]
/// for the corrections \f$ \underline g(y) \f$ to the current approximation
/// to the solution \f$ {\underline F}(y) \f$. The numerical scheme can be
/// run for any given distribution of nodes and can adapt the nodal
/// positions based on residual evaluations (including both refinement and unrefinement).

#ifndef ODE_BVP_H
#define ODE_BVP_H

#include <DenseVector.h>
#include <DenseMatrix.h>
#include <Equation_1matrix.h>
#include <OneD_Node_Mesh.h>
#include <Uncopyable.h>
#include <Residual.h>
#include <Timer.h>
#include <ArcLength_base.h>
#include <BandedLinearSystem.h>
#include <Utility.h>
#include <LinearEigenSystem_base.h>
#include <DenseLinearEigenSystem.h>

namespace CppNoddy {

  /// A templated object for real/complex vector system
  /// of first-order ordinary differential equations. The
  /// class is double templated, with the first type associated
  /// with real/complex data, and second (real/complex) type associated with
  /// a problem on the real line or line in the complex plane.
  template < typename _Type, typename _Xtype = double >
  class ODE_BVP : public ArcLength_base<_Type> {
   public:

    // we want access to the base class init_arc method even though it
    // is hidden by the derived class' init_arc.
    using ArcLength_base<_Type>::init_arc;

    /// The class is defined by a vector function for the system.
    /// \param ptr_to_equation A pointer to an Equation_1matrix object.
    /// \param nodes A vector that defines the nodal positions.
    /// \param ptr_to_left_residual A pointer to a residual object that defines the LHS boundary conditions.
    /// \param ptr_to_right_residual A pointer to a residual object that defines the RHS boundary conditions.
    ODE_BVP(Equation_1matrix<_Type, _Xtype >* ptr_to_equation,
            const DenseVector<_Xtype> &nodes,
            Residual<_Type>* ptr_to_left_residual,
            Residual<_Type>* ptr_to_right_residual);

    /// Destructor
    virtual ~ODE_BVP();

    /// A virtual method that is called prior to the linear solve
    /// stage of the solve2() method. This is called prior to any
    /// linear solve to allow the user to manually tune the matrix
    /// problem directly.
    /// \param a The Jacbian matrix to be passed to the linear solver
    /// \param b The residual vector to be passed to the linear solver
    virtual void actions_before_linear_solve(BandedMatrix<_Type>& a, DenseVector<_Type>& b)
    {}

    /// Formulate and solve the ODE using Newton iteration and
    /// a second-order finite difference scheme. The solution is
    /// stored in the publicly accessible 'solution' member data.
    void solve2();

    /// Adapt the computational mesh ONCE. We step through each interior
    /// node and evaluate the residual at that node. If the residual
    /// is less than the convergence tolerance, the node is removed.
    /// If the residual is greater than the adapt_tol parameter, then
    /// the mesh is adapted with an addition node place either side
    /// of the evaluation node.
    /// \param adapt_tol The residual tolerance at a nodal point that
    /// will lead to the mesh adaptation.
    /// \return A pair of values indicating the number of refinements
    /// and unrefinements.
    std::pair< unsigned, unsigned > adapt(const double& adapt_tol);

    /// Adaptively solve the system until no refinements or unrefinements
    /// are applied.
    /// \param adapt_tol The residual tolerance at a nodal point that
    /// will lead to the mesh adaptation.
    void adapt_until(const double& adapt_tol);

    /// Set the flag that determines if the determinant will be monitored
    /// The default is to monitor.
    /// \param flag The boolean value that the flag will be set to
    void set_monitor_det(bool flag);

    /// Initialise the class ready for arc length continuation. The base
    /// class requires a vector, so we wrap the base class method here
    /// so that the vector can be extracted from the mesh member data.
    /// \param p The pointer to the parameter
    /// \param length The initial arc length step to be taken (all in the
    ///     parameter.
    /// \param max_length The maximum arc length step to be allowed.
    void init_arc(_Type* p, const double& length, const double& max_length);

    /// Arc-length solve the system. Before this can be called the
    /// arc_init method should have been called in order to ensure
    /// we know a solution and have derivatives w.r.t. the arc-length
    /// parameter.
    double arclength_solve(const double& step);

    /// \return A handle to the solution mesh
    OneD_Node_Mesh<_Type, _Xtype>& solution();

    /// Access method to the tolerance
    /// \return A handle to the private member data TOLERANCE
    double& tolerance() {
      return TOL;
    }

    /// Access method to the maximum number of iterations
    /// \return A handle to the private member data MAX_ITERATIONS
    int& max_itns() {
      return MAX_ITERATIONS;
    }

   private:

    /// Solve the system for an initial guess by Newton iteration, this
    /// method is inherited from the ArcLength_base class and this is a
    /// wrapper that points to the 2nd-order finite-difference solution
    /// method. Because the arclength base class is for general vector
    /// systems, we pre/post convert the solution mesh to a vector.
    /// \param state An initial guess vector, the solution overwrites it.
    void solve(DenseVector<_Type>& state);

    /// Assemble the Jacobian matrix and residual vector for the
    /// BVP using the equation object and boundary condition residuals
    /// for this instance. This routine is separated into a separate
    /// method because it is used in both the solve2() and the
    /// arclength_solve() routines.
    /// \param a The matrix to be filled with the banded Jacobian. It must
    /// be sized appropriately on entry.
    /// \param b The vector to be filled with the residuals at each nodal
    /// point in the mesh.
    void assemble_matrix_problem(BandedMatrix<_Type>& a, DenseVector<_Type>& b);

    /// The solution mesh
    OneD_Node_Mesh<_Type, _Xtype> SOLUTION;
    /// maximum number of iterations to be taken
    int MAX_ITERATIONS;
    /// tolerance for convergence
    double TOL;
    /// The sign of the determinant of the Jacobian matrix for the last
    /// converged solution.
    int LAST_DET_SIGN;
    /// The equation object associated with this instance.
    Equation_1matrix<_Type, _Xtype > *p_EQUATION;
    /// Pointer to the residual defining the LHS BC
    Residual<_Type > *p_LEFT_RESIDUAL;
    /// Pointer to the residual defining the RHS BC
    Residual<_Type > *p_RIGHT_RESIDUAL;
    /// boolean flag to decide if we should monitor the determinant
    bool MONITOR_DET;

#ifdef TIME
    /// Timers for debug use
    Timer T_ASSEMBLE;
    Timer T_SOLVE;
    Timer T_REFINE;
#endif

  }
  ; // end class

  template <typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::set_monitor_det(bool flag) {
    MONITOR_DET = flag;
  }

  template <typename _Type, typename _Xtype>
  void ODE_BVP<_Type, _Xtype>::init_arc(_Type* p,
                                        const double& length,
                                        const double& max_length) {
    DenseVector<_Type> state(SOLUTION.vars_as_vector());
    this -> init_arc(state, p, length, max_length);
  }

  template <typename _Type, typename _Xtype>
  inline OneD_Node_Mesh<_Type, _Xtype>& ODE_BVP<_Type, _Xtype>::solution() {
    return SOLUTION;
  }

} // end namespace

#endif
