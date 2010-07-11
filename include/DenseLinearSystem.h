/// \file DenseLinearSystem.h
/// Specification of the linear system class.

#ifndef DENSELINEARSYSTEM_H
#define DENSELINEARSYSTEM_H

#include <DenseMatrix.h>
#include <DenseVector.h>
#include <LinearSystem_base.h>

namespace CppNoddy
{

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for dense typed problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class DenseLinearSystem : public LinearSystem_base
  {

  protected:

    typedef typename DenseMatrix<_Type>::row_iter row_iter;
    typedef typename DenseMatrix<_Type>::row_riter row_riter;
    typedef typename DenseMatrix<_Type>::row_citer row_citer;

    typedef typename DenseMatrix<_Type>::elt_iter elt_iter;
    typedef typename DenseMatrix<_Type>::elt_riter elt_riter;
    typedef typename DenseMatrix<_Type>::elt_citer elt_citer;

  public:

    /// Constructor for a dense linear system object.
    /// \param p_A A pointer to the 'A matrix', an nxn double/complex dense matrix
    /// \param p_B A pointer to the 'B vector' a size n double/complex dense vector
    /// \param which A string that indicates which solver to use
    DenseLinearSystem( DenseMatrix<_Type>* p_A, DenseVector<_Type>* p_B, std::string which = "native" );

    /// Destructor for a linear system object.
    ~DenseLinearSystem()
    {}

    /// Solve the sparse system
    void solve();

  private:

    /// Solve the linear system using LAPACK's LU solver
    void solve_lapack();

    /// Solve the linear system using the native elimination
    void solve_native();

    /// A wrapped up less_than check that throws an exception
    /// if the matrix elements are complex.
    bool lt( _Type value ) const;

    /// Back substitution routine for dense systems.
    void backsub() const;

    /// Compute the signature of the permutation vector
    /// This is used to get the sign of the det after the
    /// LAPACK routine has been called.
    /// \param pivots The pivot permutation vector
    /// \return The sign of the permutation
    int signature( const std::vector<int> &pivots ) const;

    /// A lower bound on pivot size for the native elimination routine
    double MIN_PIV;

    /// pointers to the associated containers
    DenseMatrix<_Type> *p_A;
    DenseVector<_Type> *p_B;

  };

  template <>
  inline bool DenseLinearSystem<double>::lt( double value ) const
  {
    return ( value < 0 );
  }

} //end namepsace
#endif
