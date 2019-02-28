/// \file DenseLinearEigenSystem.h
/// Specification of the dense linear eigensystem class.
/// This class links to LAPACK to perform the solver phase.

#ifndef DENSELINEAREIGENSYSTEM_BASE_H
#define DENSELINEAREIGENSYSTEM_BASE_H

#include <Types.h>
#include <LinearEigenSystem_base.h>

namespace CppNoddy {

  /// A linear Nth-order generalised eigensystem class.
  /// Here we can construct a linear eigenproblem in the form
  /// \f[ A_{NxN} \,{\underline x}_i = \lambda_i\, B_{NxN}\, {\underline x}_i \f]
  /// for dense double/complex matrices \f$ A \f$ and \f$ B \f$. The eigenvalues
  /// and eigenvectors can be tagged and retrieved as required.
  template <typename _Type>
  class DenseLinearEigenSystem : public LinearEigenSystem_base {

   public:

    /// Constructor for a linear system object.
    /// \param Aptr A pointer to a typed A matrix
    /// \param Bptr A pointer to a typed B matrix
    DenseLinearEigenSystem(DenseMatrix<_Type>* Aptr, DenseMatrix<_Type>* Bptr);

    /// Destructor for a linear system object.
    ~DenseLinearEigenSystem();

    /// Solve the matrix linear eigensystem
    void eigensolve();

    // BECAUSE OF THE WAY THE DENSE LAPACK QZ ROUTINE
    // RETURNS EIGENVALUES (in alpha/beta format) WE WILL NOT USE THE BASE
    // CLASS'S TAGGING METHODS, BUT INSTEAD DEFINE NEW VERSIONS.

    /// Tag those eigenvalues that are to the right of a specified
    /// point.
    /// \param val Tags are added or removed for val +ve or -ve
    void tag_eigenvalues_right(const int &val);

    /// Tag those eigenvalues that are to the left of a specified
    /// point.
    /// \param val Tags are added or removed for val +ve or -ve
    void tag_eigenvalues_left(const int &val);

    /// Tag those eigenvalues that are in the upper half-plane above a
    /// specified point.
    /// \param val Tags are added or removed for val +ve or -ve
    void tag_eigenvalues_upper(const int &val);

    /// Tag those eigenvalues that are in the lower half-plane below a specified
    /// point.
    /// \param val Tags are added or removed for val +ve or -ve
    void tag_eigenvalues_lower(const int &val);

    /// Tag those eigenvalues that are within a disc centred
    /// at a point in the complex plane.
    /// \param val Tags are added or removed for val +ve or -ve
    /// \param radius The radius of the disc to be considered
    void tag_eigenvalues_disc(const int &val,
                              const double &radius);

    /// Get the the tagged eigenvalues. All of the tagged eigenvalues
    /// are returned in a complex vector, with no ordering guaranteed.
    /// \return The tagged eigenvalues as a complex vector
    DenseVector<D_complex> get_tagged_eigenvalues() const;

    /// Get the the tagged eigenvectors. All of the eigenvectors associated
    /// with the tagged eigenvalues are returned, with the i-th eigenvector corresponding
    /// to the i-th eigenvalue as returned by the get_tagged_eigenvalues method.
    /// The i-th eigenvector is returned in row i of the complex dense matrix.
    /// \return The tagged eigenvectors as a complex matrix
    DenseMatrix<D_complex> get_tagged_eigenvectors() const;

   private:

    /// Solve the generalised eigenproblem and compute eigenvectors
    void eigensolve_lapack_with_vectors();

    /// Solve the generalised eigenproblem without eigenvectors
    void eigensolve_lapack_without_vectors();

    /// storage for eigenvectors and eigenvalues
    DenseVector<D_complex> m_eigenvalues_alpha;
    DenseVector<D_complex> m_eigenvalues_beta;

    /// pointers to the associated matrices
    DenseMatrix<_Type>* m_pA;
    DenseMatrix<_Type>* m_pB;

  };

} //end namepsace
#endif
