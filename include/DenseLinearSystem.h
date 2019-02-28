/// \file DenseLinearSystem.h
/// Specification of the linear system class.

#ifndef DENSELINEARSYSTEM_H
#define DENSELINEARSYSTEM_H

#include <DenseMatrix.h>
#include <DenseVector.h>

namespace CppNoddy {

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for dense typed problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class DenseLinearSystem  {

   protected:

    typedef typename DenseMatrix<_Type>::row_iter row_iter;
    typedef typename DenseMatrix<_Type>::row_riter row_riter;
    typedef typename DenseMatrix<_Type>::row_citer row_citer;

    typedef typename DenseMatrix<_Type>::elt_iter elt_iter;
    typedef typename DenseMatrix<_Type>::elt_riter elt_riter;
    typedef typename DenseMatrix<_Type>::elt_citer elt_citer;

   public:

    /// Constructor for a dense linear system object.
    /// \param m_pA A pointer to the 'A matrix', an nxn double/complex dense matrix
    /// \param m_pB A pointer to the 'B vector' a size n double/complex dense vector
    /// \param which A string that indicates which solver to use
    DenseLinearSystem(DenseMatrix<_Type>* m_pA, DenseVector<_Type>* m_pB, std::string which = "native");

    /// Destructor for a linear system object.
    ~DenseLinearSystem()
    {}

    /// Solve the sparse system
    void solve();

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

    /// Solve the linear system using the native elimination
    void solve_native();

    /// A wrapped up less_than check that throws an exception
    /// if the matrix elements are complex.
    bool lt(_Type value) const;

    /// Back substitution routine for dense systems.
    void backsub() const;

    /// Compute the signature of the permutation vector
    /// This is used to get the sign of the det after the
    /// LAPACK routine has been called.
    /// \param pivots The pivot permutation vector
    /// \return The sign of the permutation
    int signature(const std::vector<int> &pivots) const;

    /// A lower bound on pivot size for the native elimination routine
    double m_minPivot;

    /// pointers to the associated containers
    DenseMatrix<_Type> *m_pA;
    DenseVector<_Type> *m_pB;

    /// a string ID to pick out the appropriate solver
    std::string m_version;

    /// the sign of the determinant of the last solved system LHS
    int m_detSign;

    /// a flag that determines of the determinant sign should be monitored
    bool m_monitorDet;

    
  };

  template <>
  inline bool DenseLinearSystem<double>::lt(double value) const {
    return (value < 0);
  }

} //end namepsace
#endif
