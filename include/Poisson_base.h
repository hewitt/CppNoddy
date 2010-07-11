/// \file Poisson_base.h
/// Specification of a Poisson base class ... assuming Dirichlet
/// conditions on the boundary.

#ifndef POISSON_BASE_H
#define POISSON_BASE_H

#include <DenseMatrix.h>
#include <BandedMatrix.h>
#include <Uncopyable.h>
#include <BandedLinearSystem.h>

namespace CppNoddy
{

  /// An base class for simple Poisson problems

  class Poisson_base : private Uncopyable
  {

  public:

    /// A simple Finite Differenced Poisson problem class.
    /// \param left The position of the left-hand boundary.
    /// \param right The position of the right-hand boundary.
    /// \param bottom The position of the bottom boundary.
    /// \param top The position of the top boundary.
    /// \param Nx Number of nodes in "x-direction".
    /// \param Ny Number of nodes in "y-direction".
    /// \param source_ptr A pointer to the (Nx by Ny) matrix that is to be the
    ///     source for the Poisson problem.
    Poisson_base( const double& left,
                  const double& right,
                  const double& bottom,
                  const double& top,
                  const unsigned& Nx,
                  const unsigned& Ny,
                  DenseMatrix<double>* source_ptr )
    {
      this -> NX = Nx;
      this -> NY = Ny;
      this -> DX = ( right - left ) / ( Nx - 1 );
      this -> DY = ( top - bottom ) / ( Ny - 1 );
      this -> LEFT = left;
      this -> RIGHT = right;
      this -> TOP = top;
      this -> BOTTOM = bottom;
      this -> p_SOURCE = source_ptr;
      A = BandedMatrix<double>( Nx * Ny, Nx, 0.0 );
      B = DenseVector<double>( Nx * Ny, 0.0 );
#ifdef DEBUG
      std::cout << "[DEBUG] Poisson problem with a " << NX << " x " << NY << " mesh\n";
      std::cout << "[DEBUG] DX = " << DX << " DY = " << DY << "\n";
      std::cout << "[DEBUG] domain is [" << LEFT << ", " << RIGHT << "] x [" << BOTTOM << ", " << TOP << "]\n";
#endif
    }

    virtual ~Poisson_base()
    {
    }

    /// Direct solver for the whole 2D domain.
    virtual void solve()
    {
      std::string problem;
      problem = "The Poisson_base::solve method has not been implemented.\n";
      problem += "You have to implement this method to define the Poisson solver.\n";
      throw ExceptionRuntime( problem );
    }

    /// A dump function for data output.
    virtual void dump()
    {
      for ( unsigned i = 0; i < NX; ++i )
      {
        for ( unsigned j = 0; j < NY; ++j )
        {
          std::cout << LEFT + i * DX << " "
                    << BOTTOM + j * DY << " "
                    << p_SOURCE -> get( i, j ) << "\n";
        }
        std::cout << "\n";
      }
    }

    /// Assemble the LHS of the linear Poisson matrix problem
    virtual void assemble_LHS()
    {
      std::string problem;
      problem = "The Poisson_base::assemble_LHS method has not been implemented.\n";
      problem += "You have to implement this method to define the Poisson solver.\n";
      throw ExceptionRuntime( problem );
    }

  protected:

    /// geometry of the domain/meesh
    unsigned NX, NY;
    double DX, DY;
    double LEFT, RIGHT, TOP, BOTTOM;
    /// pointer to the source matrix
    DenseMatrix<double>* p_SOURCE;
    /// The containers for the banded matrix problem
    BandedMatrix<double> A;
    DenseVector<double> B;
    /// The banded linear system
    BandedLinearSystem<double>* p_SYSTEM;

  };

}

#endif // POISSON_BASE_H
