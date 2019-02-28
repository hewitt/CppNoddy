/// \file reversed_BL.h
/// An implementation of the zig-zag modification for unsteady parabolic
/// marching with reverse flow. The reversal detection is hard-wired (hacked)
/// in the method below -- applicable to a box scheme boundary layer approach.

#ifndef REVERSED_BL_H
#define REVERSED_BL_H

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
  class reversed_BL : private Uncopyable {
   public:

    /// The class is defined by a vector function for the system.
    /// \param equation_ptr A pointer to an inherited Equation object.
    /// \param xnodes A vector that defines the nodal x-positions.
    /// \param ynodes A vector that defines the nodal y-positions.
    /// \param ptr_to_bottom_residual A pointer to a residual object that defines the y=y1 boundary conditions.
    /// \param ptr_to_top_residual A pointer to a residual object that defines the y=y2 boundary conditions.
    reversed_BL(Equation_3matrix<_Type > *equation_ptr,
                const DenseVector<double>& xnodes,
                const DenseVector<double>& ynodes,
                Residual_with_coords<_Type>* ptr_to_bottom_residual,
                Residual_with_coords<_Type>* ptr_to_top_residual);

    /// Destructor
    ~reversed_BL();

    /// Copy the current solution to the previous solution
    void update_previous_solution() {
      PREV_SOLN = SOLN;
      UPDATED = true;
    }

    /// A Crank-Nicolson 'time' stepper.
    void step2(const double& dt);

    /// A Crank-Nicolson 'time' stepper.
    void bidirectional_step2(const double& dt);

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
    /// \param j The index of the last x-position for which a solution is known, hence the matrix
    /// assembled is for the x=x_{j+1} position.
    void assemble_matrix_problem(BandedMatrix<_Type>& a, DenseVector<_Type>& b, const std::size_t& j);


    void first_order_assemble_matrix_problem(BandedMatrix<_Type>& a, DenseVector<_Type>& b, const std::size_t& j);


    void bidirectional_assemble_matrix_problem(BandedMatrix<_Type>& a, DenseVector<_Type>& b, const std::size_t& j);

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
  inline TwoD_Node_Mesh<_Type>& reversed_BL<_Type>::solution() {
    return SOLN;
  }

  template <typename _Type>
  inline double& reversed_BL<_Type>::t() {
    return T;
  }

} // end namespace

#endif
