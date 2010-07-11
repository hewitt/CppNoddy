/// \file BandedLinearEigenSystem.h
/// Specification of the banded linear eigensystem class.
/// This class links to ARPACK to perform the solver phase.

#ifndef BANDEDLINEAREIGENSYSTEM_BASE_H
#define BANDEDLINEAREIGENSYSTEM_BASE_H

#include <set>

#include <Uncopyable.h>
#include <Types.h>
#include <LinearEigenSystem_base.h>
#include <Utility.h>

namespace CppNoddy
{

  /// A linear Nth-order generalised eigensystem class.
  /// Here we can construct a linear eigenproblem in the form
  /// \f[ A_{NxN} \,{\underline x}_i = \lambda_i\, B_{NxN}\, {\underline x}_i \f]
  /// for Banded double/complex matrices \f$ A \f$ and \f$ B \f$. The eigenvalues
  /// and eigenvectors can be tagged and retrieved as required.
  template <typename _Type>
  class BandedLinearEigenSystem : public LinearEigenSystem_base
  {

  public:

    /// Constructor for a linear system object.
    /// \param Aptr A pointer to a typed A matrix
    /// \param Bptr A pointer to a typed B matrix
    BandedLinearEigenSystem( BandedMatrix<_Type>* Aptr, BandedMatrix<_Type>* Bptr );

    /// Destructor for a linear system object.
    ~BandedLinearEigenSystem();

    /// Access the number of Arnoldi vectors.
    /// \return A hande to the number of Arnoldi vectors
    unsigned& access_narnoldi();

    /// Access the number of eigenvalues.
    /// \return A hande to the number of eigenvalues
    unsigned& access_nev();

    /// Solve the matrix linear eigensystem
    void eigensolve();


  private:

    /// Solve the generalised eigenproblem and compute eigenvectors
    void eigensolve_arpack();

    /// pointer to the LHS matrix
    BandedMatrix<_Type>* p_A;
    /// pointer to the RHS matrix
    BandedMatrix<_Type>* p_B;

    /// the number of Arnoldi vectors to use
    unsigned NARNOLDI;

    /// the target number of eigenvalues to get
    unsigned NEV;

  };

} //end namepsace
#endif
