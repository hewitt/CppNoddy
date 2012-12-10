/// \file Poisson_Cartesian.cpp
/// Implementation of the Poisson Cartesian class.

#include <string>

#include <Types.h>
#include <Poisson_Cartesian.h>
#include <Poisson_base.h>
#include <Exceptions.h>
#include <BandedLinearSystem.h>
#include <SparseLinearSystem.h>
#include <Utility.h>

namespace CppNoddy
{


  Poisson_Cartesian::Poisson_Cartesian( const double& left,
                                        const double& right,
                                        const double& bottom,
                                        const double& top,
                                        const unsigned& Nx,
                                        const unsigned& Ny,
                                        DenseMatrix<double>* source_ptr ) :
      Poisson_base( left, right, bottom, top, Nx, Ny, source_ptr )
  {
#ifdef LAPACK
    assemble_LHS();
    p_SYSTEM = new BandedLinearSystem<double>( &A, &B, "lapack" );

    // since the LHS of the problem is fixed, we can LU decompose
    // here in the constructor for a zero RHS, then any actual solve is
    // done just done as a re_solve using the decomposition
    // this is obviously much faster! This requires the LAPACK solver
    // since the native solver does not do LU.
    p_SYSTEM -> solve();
#else
    p_SYSTEM = new BandedLinearSystem<double>( &A, &B, "native" );
#endif
  }

  void Poisson_Cartesian::solve()
  {
    // top & bottom
    for ( unsigned i = 0; i < NX; ++i )
    {
      B[ i ] = p_SOURCE -> get( i, 0 );
      B[ ( NY - 1 ) * NX + i ] = p_SOURCE -> get( i, NY - 1 );
    }

    // left & right
    for ( unsigned j = 0; j < NY; ++j )
    {
      B[ j * NX ] = p_SOURCE -> get( 0, j );
      B[ j * NX + NX - 1 ] = p_SOURCE -> get( NX - 1, j );
    }

    // interior points
    for ( unsigned i = 1; i < NX - 1; ++i )
    {
      for ( unsigned j = 1; j < NY - 1; ++j )
      {
        B[ j * NX + i ] = p_SOURCE -> get( i, j );
      }
    }
#ifdef LAPACK
    // if we have LAPACK, use the existing LU decomposition
    p_SYSTEM -> re_solve_lapack();
#else
    // native solver sets up the whole linear system from scratch & solves
    assemble_LHS();
    p_SYSTEM -> solve();
#endif

    for ( unsigned i = 0; i < NX; ++i )
    {
      for ( unsigned j = 0; j < NY; ++j )
      {
        p_SOURCE -> set( i, j ) = B[ j * NX + i ];
      }
    }
  }

  void Poisson_Cartesian::assemble_LHS()
  {
    const double inv_dx2 = 1 / ( DX * DX );
    const double inv_dy2 = 1 / ( DY * DY );
    // clear the matrix
    Utility::fill( A, 0.0 );
    // top & bottom Dirichlet
    for ( unsigned i = 0; i < NX; ++i )
    {
      A( i, i ) = 1.0;
      A( ( NY - 1 ) * NX + i, ( NY - 1 ) * NX + i ) = 1.0;
    }
    // left & right Dirichlet
    for ( unsigned j = 0; j < NY; ++j )
    {
      A( j * NX, j * NX ) = 1.0;
      A( j * NX + NX - 1, j * NX + NX - 1 ) = 1.0;
    }
    // interior points
    for ( unsigned i = 1; i < NX - 1; ++i )
    {
      for ( unsigned j = 1; j < NY - 1; ++j )
      {
        A( j * NX + i, j * NX + i ) = -2 * ( inv_dx2 + inv_dy2 );
        A( j * NX + i, j * NX + i - 1 ) = inv_dx2;
        A( j * NX + i, j * NX + i + 1 ) = inv_dx2;
        A( j * NX + i, j * NX + i - NX ) = inv_dy2;
        A( j * NX + i, j * NX + i + NX ) = inv_dy2;
      }
    }
  }


}   // end namespace CppNoddy
