/// \file Poisson_Cartesian.h
/// Specification of a Poisson problem in a 2-D
/// Cartesian coordinate system with Dirichlet
/// boundary conditions -- you *really* don't want
/// to use this unless you link to LAPACK.

#ifndef POISSON_CARTESIAN_H
#define POISSON_CARTESIAN_H

#include <DenseMatrix.h>
#include <Poisson_base.h>
#include <BandedLinearSystem.h>

namespace CppNoddy
{

  /// An object for Cartesian Poisson problems with Dirichlet boundary
  /// conditions
  class Poisson_Cartesian : public Poisson_base
  {

  public:

    /// Constructor for a simple Finite Differenced Poisson problem class.
    /// \param left The x-position of the left-hand boundary.
    /// \param right The x-position of the right-hand boundary.
    /// \param bottom The y-position of the bottom boundary.
    /// \param top The y-position of the top boundary.
    /// \param Nx Number of nodes in x-direction.
    /// \param Ny Number of nodes in y-direction.
    /// \param source_ptr A pointer to a DenseMatrix<double> source term.
    Poisson_Cartesian( const double& left,
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

  };

}

#endif // POISSON_CARTESIAN_H
