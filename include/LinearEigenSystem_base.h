/// \file LinearEigenSystem_base.h
/// Specification of the linear eigensystem base class.

#ifndef LINEAREIGENSYSTEM_BASE_H
#define LINEAREIGENSYSTEM_BASE_H

#include <set>

#include <Types.h>

namespace CppNoddy {

  /// A linear Nth-order generalised eigensystem base class.
  /// Here we can construct a linear eigenproblem in the form
  /// \f[ A_{NxN} \,{\underline x}_i = \lambda_i\, B_{NxN}\, {\underline x}_i \f]
  /// for matrices \f$ A \f$ and \f$ B \f$. The eigenvalues
  /// and eigenvectors can be tagged and retrieved as required.
  class LinearEigenSystem_base {

   protected:

    typedef std::set
    < unsigned, std::less<unsigned> >::iterator iter;
    typedef std::set
    < unsigned, std::less<unsigned> >::const_iterator citer;

   public:

    /// Constructor for a linear system object.
    LinearEigenSystem_base();

    /// Destructor for a linear system object.
    virtual ~LinearEigenSystem_base();

    /// Solve the matrix linear eigensystem
    virtual void eigensolve();

    /// Compute the eigenvectors in any eigenvalue computation
    /// \param flag The boolean value to set.
    void set_calc_eigenvectors(bool flag);

    /// Set the shift value to be used in tagging
    /// \param z The shift value to be used
    void set_shift(const D_complex& z);

    /// Get the shift value associated with this class (used for tagging)
    /// and eigenvalue tagging
    D_complex get_shift() const;

    /// Tag those eigenvalues that are to the right of a specified
    /// shft point.
    /// \param val Tags are added or removed for val +ve or -ve
    virtual void tag_eigenvalues_right(const int &val);

    /// Tag those eigenvalues that are to the left of a specified
    /// point.
    /// \param val Tags are added or removed for val +ve or -ve
    virtual void tag_eigenvalues_left(const int &val);

    /// Tag those eigenvalues that are in the upper half-plane above a
    /// specified point.
    /// \param val Tags are added or removed for val +ve or -ve
    virtual void tag_eigenvalues_upper(const int &val);

    /// Tag those eigenvalues that are in the lower half-plane below a specified
    /// point.
    /// \param val Tags are added or removed for val +ve or -ve
    virtual void tag_eigenvalues_lower(const int &val);

    /// Tag those eigenvalues that are within a disc centred
    /// at a point in the complex plane.
    /// \param val Tags are added or removed for val +ve or -ve
    /// \param radius The radius of the disc to be considered
    virtual void tag_eigenvalues_disc(const int &val, const double &radius);

    /// Tag all the eigenvalues returned from the ARPACK routine.
    virtual void tag_eigenvalues_all();

    virtual void clear_all_tags() {
      m_tagged_indices.clear();
    }

    /// Get the the tagged eigenvalues. All of the tagged eigenvalues
    /// are returned in a complex vector, with no ordering guaranteed.
    /// \return The tagged eigenvalues as a complex vector
    virtual DenseVector<D_complex> get_tagged_eigenvalues() const;

    /// Get the the tagged eigenvectors. All of the eigenvectors associated
    /// with the tagged eigenvalues are returned, with the i-th eigenvector corresponding
    /// to the i-th eigenvalue as returned by the get_tagged_eigenvalues method.
    /// The i-th eigenvector is returned in row i of the complex dense matrix.
    /// \return The tagged eigenvectors as a complex matrix
    virtual DenseMatrix<D_complex> get_tagged_eigenvectors() const;

   protected:

    /// storage for eigenvectors and eigenvalues
    DenseVector<D_complex> m_all_eigenvalues;
    DenseMatrix<D_complex> m_all_eigenvectors;

    /// a set of tagged eigenvalue indices
    std::set
    < unsigned, std::less<unsigned> > m_tagged_indices;

    /// The complex shift value
    D_complex m_shift;

    /// calculate the eigenvectors in any eigenvalue computation
    bool m_calc_eigenvectors;

  };

} //end namepsace
#endif
