/// \file SparseLinearEigenSystem.h
/// Specification of the sparse linear eigensystem class.
/// This class links to complex SLEPc to perform the solver phase.

#ifndef SPARSELINEAREIGENSYSTEM_BASE_H
#define SPARSELINEAREIGENSYSTEM_BASE_H

#include <set>
#include <SparseMatrix.h>
#include <Uncopyable.h>
#include <Types.h>
#include <LinearEigenSystem_base.h>

#ifdef SLEPC
#include <slepc.h>

namespace CppNoddy {
  PetscErrorCode monitor_function(EPS eps, PetscInt its,PetscInt nconv,
                                  PetscScalar *eigr,PetscScalar *eigi,
                                  PetscReal* errest,
                                  PetscInt nest,void *mctx) {
    std::cout << "[MONITOR] nconv = " << nconv;
    std::cout << " its = " << its << "\n";
    std::cout << "[MONITOR] est_err() = ";
    for(int i = 0; i < nest; i++) {
      std::cout << errest[i] << " ";
    }
    std::cout << "\n";
    return 0;
  }

  /// A linear Nth-order generalised eigensystem class.
  /// Here we can construct a linear eigenproblem in the form
  /// \f[ A_{NxN} \,{\underline x}_i = \lambda_i\, B_{NxN}\, {\underline x}_i \f]
  /// for Banded double/complex matrices \f$ A \f$ and \f$ B \f$. The eigenvalues
  /// and eigenvectors can be tagged and retrieved as required.
  template <typename _Type>
  class SparseLinearEigenSystem : public LinearEigenSystem_base {

   public:

    /// Constructor for a linear system object.
    /// \param Aptr A pointer to a typed A matrix
    /// \param Bptr A pointer to a typed B matrix
    SparseLinearEigenSystem(SparseMatrix<_Type>* Aptr, SparseMatrix<_Type>* Bptr);

    /// Destructor for a linear system object.
    ~SparseLinearEigenSystem();

    /// Access the (actual) number of (converged) eigenvalues found.
    /// \return A hande to the number of eigenvalues
    unsigned get_nconv() const;

    /// Request a certain number of eigenvalues (not guaranteed)
    /// \param n The number of eigenvalues being asked for
    void set_nev(unsigned n);

    /// Set target for the shift-invert algorithm.
    /// \param target The target eigenvalue
    void set_target(std::complex<double> target);

    /// Defines the ordering of returned set of eigenvalues/vectors
    /// \param order_string A string that defines the SLEPc enum options
    ///  for ordering (see EPSWhich)
    void set_order(std::string order_string);

    /// Gives a handle to the boolean that sets a region in the complex plane
    bool& region_defined();

    /// Set a rectangular region of the complex plane in which to look for eigenvalues.
    /// Calling this will also set REGION_DEFINED to true.
    /// \param a Real=a defines left edge of the rectangular region
    /// \param b Real=b defines right edge of the rectangular region
    /// \param c Imag=c defines bottom edge of the rectangular region
    /// \param d Imag=d defines top edge of the rectangular region
    void set_region(const double& a, const double& b, const double& c, const double& d);

    /// Gives a handle to the boolean that sets if an initial guess has been used
    bool& guess_defined();

    /// Set the initial guess
    /// \param guess The vector of the guess
    void set_initial_guess(const DenseVector<_Type>& guess);

    /// Solve the matrix linear eigensystem
    void eigensolve();

   private:

    /// Solve the generalised eigenproblem and compute eigenvectors
    void eigensolve_slepc();

    /// number of eigenvalues requested
    unsigned m_nev;
    /// number of (converged) eigenvalues located by the solver will be changed by the solver.
    unsigned m_nconv;

    /// defines how to seek and order eigenvalues (e.g. smallest real part?)
    EPSWhich m_order;

    /// defines if ev's are to be sought in a specific region
    bool m_region_defined;
    /// defines if an initial guess is to be used
    bool m_guess_defined;

    /// four numbers that define a rectangular region in the complex plane
    double m_real_left,m_real_right,m_imag_bottom,m_imag_top;

    /// stores an initial guess to work with
    DenseVector<_Type> m_initial_guess;

    /// pointer to the LHS matrix
    SparseMatrix<_Type>* m_pA;
    /// pointer to the RHS matrix
    SparseMatrix<_Type>* m_pB;

    // SLEPc* p_LIBRARY;
  };

} //end namepsace
#endif


#endif
