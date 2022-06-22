/// \file BandedLinearSystem.h
/// Specification of the linear system class.

#ifndef BANDEDLINEARSYSTEM_H
#define BANDEDLINEARSYSTEM_H

#include <Types.h>

namespace CppNoddy {

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for dense typed problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class BandedLinearSystem {

   public:

    /// Constructor for a banded linear system object.
    /// \param Aptr A pointer to the 'A matrix', an NxN double/complex banded matrix
    /// \param Bptr A pointer to the 'B vector' a size N double/complex dense vector
    /// \param which A string that indicates which solver to use
    BandedLinearSystem(BandedMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which = "native");

    /// Destructor for a linear system object.
    ~BandedLinearSystem()
    {}

    /// Solve the banded system
    void solve();

    /// Resolve the banded system
    void re_solve_lapack();

    /// Get the sign of the determinant of the LHS matrix
    /// from the linear system just computed.
    /// \return The sign of the determinant of the
    /// LAST solved system.
    int get_det_sign() const;

    /// Store the sign of the determinant of the LHS matrix
    /// every time a solve is requested on a real system.
    /// \param flag The boolean value to set.
    void set_monitor_det(bool flag);
    
   private:

    /// Solve the linear system using LAPACK's LU solver
    void solve_lapack();

    /// Solve the linear system using the native elimination.
    /// Note that access to the banded matrix is done through
    /// the operator() methods and this method is sloooow!
    void solve_native();

    /// A wrapped up less_than check that throws an exception
    /// if the matrix elements are complex.
    bool lt(_Type value) const;

    /// Back substitution routine for dense systems.
    /// \param A The upper triangular matrix LHS
    /// \param B The dense vector RHS
    void backsub() const;

    /// Compute the signature of the permutation vector
    int signature(const std::vector<int> &pivots) const;

    /// pointers to the associated containers
    BandedMatrix<_Type> *m_pA;
    DenseVector<_Type> *m_pB;

    /// a string ID to pick out the appropriate solver
    std::string m_version;

    /// the sign of the determinant of the last solved system LHS
    int m_detSign;

    /// a flag that determines of the determinant sign should be monitored
    bool m_monitorDet;
    
    /// to allow resolves, we store pivots
    std::vector<int> m_pivots;
  };

} //end namepsace
#endif
