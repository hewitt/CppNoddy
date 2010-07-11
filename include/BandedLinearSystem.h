/// \file BandedLinearSystem.h
/// Specification of the linear system class.

#ifndef BANDEDLINEARSYSTEM_H
#define BANDEDLINEARSYSTEM_H

#include <Types.h>
#include <LinearSystem_base.h>


namespace CppNoddy
{

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for dense typed problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  template <typename _Type>
  class BandedLinearSystem : public LinearSystem_base
  {

  public:

    /// Constructor for a banded linear system object.
    /// \param Aptr A pointer to the 'A matrix', an NxN double/complex banded matrix
    /// \param Bptr A pointer to the 'B vector' a size N double/complex dense vector
    /// \param which A string that indicates which solver to use
    BandedLinearSystem( BandedMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which = "native" );

    /// Destructor for a linear system object.
    ~BandedLinearSystem()
    {}

    /// Solve the banded system
    void solve();

    /// Resolve the banded system
    void re_solve_lapack();

  private:

    /// Solve the linear system using LAPACK's LU solver
    void solve_lapack();

    /// Solve the linear system using the native elimination.
    /// Note that access to the banded matrix is done through
    /// the operator() methods and this method is sloooow!
    void solve_native();

    /// A wrapped up less_than check that throws an exception
    /// if the matrix elements are complex.
    bool lt( _Type value ) const;

    /// Back substitution routine for dense systems.
    /// \param A The upper triangular matrix LHS
    /// \param B The dense vector RHS
    void backsub() const;

    /// Compute the signature of the permutation vector
    int signature( const std::vector<int> &pivots ) const;

    /// pointers to the associated containers
    BandedMatrix<_Type> *p_A;
    DenseVector<_Type> *p_B;

    /// to allow resolves, we store pivots
    std::vector<int> STORED_PIVOTS;
  };

} //end namepsace
#endif
