/// \file Poisson_meridional.h
/// Specification of a Poisson problem in the meridional plane
/// of a cylindrical polar coordinate system with Dirichlet
/// boundary conditions -- you *really* don't want to use
/// this unless you link to LAPACK.

#ifndef POISSON_MERIDIONAL_H
#define POISSON_MERIDIONAL_H

#include <DenseMatrix.h>
#include <Poisson_base.h>
#include <BandedLinearSystem.h>

namespace CppNoddy
{

  /// An object for Poisson problems in the meridional plane of cylindrial polars

  class Poisson_meridional : public Poisson_base
  {

  public:

    /// A simple Finite Differenced Poisson problem class.
    /// \param left The x-position of the left-hand boundary.
    /// \param right The x-position of the right-hand boundary.
    /// \param bottom The y-position of the bottom boundary.
    /// \param top The y-position of the top boundary.
    /// \param Nx Number of nodes in x-direction.
    /// \param Ny Number of nodes in y-direction.
    /// \param source_ptr A pointer to a DenseMatrix<double> source term.
    Poisson_meridional( const double& left,
                        const double& right,
                        const double& bottom,
                        const double& top,
                        const unsigned& Nx,
                        const unsigned& Ny,
                        DenseMatrix<double>* source_ptr );

    /// Direct solver for the whole 2D domain.
    void solve();

    /// Assemble the LHS of the linear Poisson matrix problem
    void assemble_LHS();

    /// For a class of Stokes streamfunction problem there is a
    /// sign change in the Laplacian-like operator. This sets
    /// the sign change.
    void set_stokes_streamfn()
    {
      STOKES_STREAMFN = -1;
      // since the operator on the LHS is different we need to re-factorise it
#ifdef LAPACK
      assemble_LHS();
      p_SYSTEM -> solve();
#endif
    }

    /// For a class of Stokes streamfunction problem there is a
    /// sign change in the Laplacian-like operator. This unsets
    /// the sign change, to give the standard Laplacian.
    void unset_stokes_streamfn()
    {
      STOKES_STREAMFN = 1;
      // since the operator on the LHS is different we need to re-factorise it
#ifdef LAPACK
      assemble_LHS();
      p_SYSTEM -> solve();
#endif
    }

  private:

    /// In some cases involving a Stokes streamfunction the operator
    /// is not standard Laplacian, but has a sign change. This value
    /// is set to -1 for such cases.
    int STOKES_STREAMFN;

  };

}

#endif // POISSON_MERIDIONAL_H
