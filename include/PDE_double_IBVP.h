/// \file PDE_double_IBVP.h
/// A specification of a class for an \f$ n^{th} \f$-order IBVP of the form
/// \f[ M_2( {\underline f}(y,x,t), y, x, t )\cdot {\underline f}_x (y,x,t) + M_1( {\underline f}(y,x,t), y,x, t )\cdot {\underline f}_t (y,x,t) + M_0( {\underline f}(y,x,t), y,x, t )\cdot {\underline f}_y (y,x,t) = {\underline R}( {\underline f}(y,x,t), y,x, t )\,, \f]
/// subject to \f$ n \f$ conditions defined at \f$ y = y_{bottom} \f$ and
/// \f$ y_{top} \f$ for some components of \f$ {\underline f}(y,x,t) \f$.
/// Here \f$ M_2 \f$ is a matrix for the x-variation and \f$ M_1 \f$ is the mass
/// matrix for the t-variation, whilst \f$ M_0 \f$ is not restricted to be an identity matrix (though it often will be).
/// The solution at the new time step \f$ t+\Delta t \f$ and spatial location \f$ x = x_j \f$ is
/// \f[ {\underline f}^{new} = {\underline F}_{j} + {\underline g} \f]
/// where \f$ {\underline F}_{j} \f$ is the current guess at the solution at this point
/// and \f$ {\underline g} \f$ is the linearised correction.
/// The solution at the previous time \f$ t \f$ at \f$ x = x_{j} \f$ is
/// \f[ {\underline f}^{old} = {\underline O}_{j} \f]
/// A Crank-Nicolson method is employed with the linearised problem at the mid-t-point
/// \f$ t + \Delta t /2 \f$ and mid-x-point \f$ ( x_{j+1} + x_{j} ) / 2 \f$ being:
/// \f[ \frac{2}{x_{j+1} - x_{j}}\, M_2 \cdot {\underline g } + \frac{2}{\Delta t}\, M_1 \cdot {\underline g } +  M_0 \cdot {\underline g}_y - J \cdot {\underline g} \f]
/// \f[ + J_{2} \cdot \frac{ ( \underline F_{j+1} + \underline O_{j+1} - \underline F_j - \underline O_j) }{ 2(x_{j+1} - x_j) } \cdot {\underline g} \f]
/// \f[ + J_{1} \cdot \frac{(\underline F_{j+1} + \underline F_j - \underline O_j - \underline O_{j+1})}{2\Delta t} \cdot {\underline g} \f]
/// \f[ + J_{0} \cdot \frac{( \underline F_{j+1} + \underline F_j + \underline O_j + \underline O_{j+1})_y}{4} \cdot {\underline g} \f]
/// \f[ = 4 {\underline R} - M_0 \cdot \biggl ({\underline F_j}_y + {\underline O_j}_y + {\underline F_{j+1}}_y + {\underline O_{j+1}}_y \biggr ) - \frac{2}{x_{j+1}-x_j} M_1 \cdot \biggl ( {\underline F}_{j+1} + {\underline O}_{j+1} - {\underline F}_{j} - {\underline O}_{j} \biggr ) - \frac{2}{\Delta t} M_2 \cdot \biggl ( {\underline F}_{j+1} + {\underline F}_{j} - {\underline O}_{j+1} - {\underline O}_{j} \biggr ) \f]
/// Where \f$ M_{0,1,2}, J, J_{0,1,2}, R \f$ are evaluated at the mid-t/mid-x step with arguments \f$ \left ( \frac{\underline F_{j} + \underline O_{j} + \underline F_{j+1} + \underline O_{j+1} }{4}, y, ( x_{j+1} + x_{j} ) / 2, t + \frac{\Delta t}{2} \right ) \f$,
/// with \f$ J_{0,1,2} \f$ denoting the Jacobian of the mass matrices \f$ \partial {M_{0,1,2}}_{ij} / \partial f_k \f$.
/// This problem is solved by second-order central differencing the equation at
/// the spatial (\f$ y \f$) inter-node mid points.

#ifndef PDE_DOUBLE_IBVP_H
#define PDE_DOUBLE_IBVP_H

#include <DenseVector.h>
#include <DenseMatrix.h>
#include <BandedMatrix.h>
#include <Equation_3matrix.h>
#include <Residual_with_coords.h>
#include <TwoD_Node_Mesh.h>
#include <Uncopyable.h>
#include <Timer.h>

namespace CppNoddy {

  /// A templated object for real/complex vector system
  /// of unsteady equations.

  template <typename _Type>
  class PDE_double_IBVP : private Uncopyable {
   public:

    /// The class is defined by a vector function for the system.
    /// \param equation_ptr A pointer to an inherited Equation object.
    /// \param xnodes A vector that defines the nodal x-positions.
    /// \param ynodes A vector that defines the nodal y-positions.
    /// \param ptr_to_bottom_residual A pointer to a residual object that defines the y=y1 boundary conditions.
    /// \param ptr_to_top_residual A pointer to a residual object that defines the y=y2 boundary conditions.
    PDE_double_IBVP(Equation_3matrix<_Type > *equation_ptr,
                    const DenseVector<double>& xnodes,
                    const DenseVector<double>& ynodes,
                    Residual_with_coords<_Type>* ptr_to_bottom_residual,
                    Residual_with_coords<_Type>* ptr_to_top_residual);

    /// Destructor
    ~PDE_double_IBVP();

    /// Copy the current solution to the previous solution
    void update_previous_solution() {
      PREV_SOLN = SOLN;
      UPDATED = true;
    }

    /// A Crank-Nicolson 'time' stepper.
    void step2(const double& dt);

    /// Return a reference to the current value of the 'timelike' coordinate
    /// \return A handle to the current timelinke coordinate stored in the object
    double& t();

    /// Return a reference to the convergence tolerance
    double& tolerance() {
      return TOL;
    }

    /// \return A handle to the solution mesh
    TwoD_Node_Mesh<_Type>& solution();

   private:

    /// Assembles the matrix problem for a BVP solve at the current time level
    /// and for the j+1 x-location.
    /// \param a The LHS (banded) matrix.
    /// \param b The RHS (dense) vector.
    /// \param j The index of the last x-position for which a solution is known.
    void assemble_matrix_problem(BandedMatrix<_Type>& a, DenseVector<_Type>& b, const std::size_t& j);

    /// The solution mesh at the current time value
    TwoD_Node_Mesh<_Type> SOLN;
    /// The solution at the previous time step
    TwoD_Node_Mesh<_Type> PREV_SOLN;
    /// tolerance
    double TOL;
    /// The current value of the timelike variable
    double T;
    /// (uniform) temporal step size
    double DT;
    /// maximum number of iterations
    int MAX_ITERATIONS;
    /// The function associated with this instance.
    Equation_3matrix<_Type > *p_EQUATION;
    /// Pointer to the residual defining the 'bottom' BC
    Residual_with_coords<_Type > *p_BOTTOM_RESIDUAL;
    /// Pointer to the residual defining the 'top' BC
    Residual_with_coords<_Type > *p_TOP_RESIDUAL;
    /// A flag to make sure the update has been done
    bool UPDATED;

#ifdef TIME
    /// Timers for debug use
    Timer T_ASSEMBLE;
    Timer T_SOLVE;
#endif

  }
  ; // end class

  template <typename _Type>
  inline TwoD_Node_Mesh<_Type>& PDE_double_IBVP<_Type>::solution() {
    return SOLN;
  }

  template <typename _Type>
  inline double& PDE_double_IBVP<_Type>::t() {
    return T;
  }

} // end namespace

#endif

