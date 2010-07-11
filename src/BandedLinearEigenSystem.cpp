/// \file BandedLinearEigenSystem.cpp
/// Implementation for the BandedLinearEigenSystem class
/// This class links to ARPACK to perform the solver phase.

#include <vector>
#include <set>
#include <algorithm>

#include <FortranARPACK.h>
#include <FortranData.h>
#include <BandedLinearEigenSystem.h>
#include <DenseLinearEigenSystem.h>
#include <Exceptions.h>
#include <Types.h>

namespace CppNoddy
{

  template <typename _Type>
  BandedLinearEigenSystem<_Type>::BandedLinearEigenSystem( BandedMatrix<_Type>* Aptr, BandedMatrix<_Type>* Bptr ) :
      LinearEigenSystem_base(),
      NEV( 4 )
  {
    p_A = Aptr;
    p_B = Bptr;
    const unsigned default_nArnoldi( 8 );
    NARNOLDI = std::min( ( unsigned ) ( p_A -> nrows() ), default_nArnoldi );
  }

  template <typename _Type>
  BandedLinearEigenSystem<_Type>::~BandedLinearEigenSystem()
  {}

  template <typename _Type>
  unsigned& BandedLinearEigenSystem<_Type>::access_narnoldi()
  {
    return NARNOLDI;
  }

  template <typename _Type>
  unsigned& BandedLinearEigenSystem<_Type>::access_nev()
  {
    return NEV;
  }

  template <typename _Type>
  void BandedLinearEigenSystem<_Type>::eigensolve()
  {
    eigensolve_arpack();
  }

  template <>
  void BandedLinearEigenSystem<double>::eigensolve_arpack()
  {
#ifndef ARPACK
    std::string problem;
    problem = "The BandedLinearEigenSystem::eigensolve_arpack_without_vectors method has been called\n";
    problem += "but the compiler option -DARPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal( problem );
#else
    // size of the system
    const unsigned N( p_A -> nrows() );
    // number of vectors to be used
    const unsigned ncv( NARNOLDI );
    if ( ( nev * 2 > NARNOLDI ) || ( nev > N ) || ( NARNOLDI > N ) )
    {
      std::string problem;
      problem = "The BandedLinearEigenSystem::eigensolve_arpack_without_vectors has been called\n";
      problem += "with inconsistent choices for the number of eigenvalues to find and/or \n";
      problem += "the NARNOLDI vectors.\n";
      throw ExceptionExternal( problem, 0 );
    }
    // leading dimension of the banded matrix
    const unsigned lda( 3 * p_A -> noffdiag() + 1 );
    // leading dimension of additional work arrays
    const unsigned ldv = N;
    const unsigned ldz = N;
    // parameter vector
    std::vector<int> iparam( 11, 0 );
    // workspaces
    DD_vector rfac( lda * N, 0.0 );
    DD_vector v( ldv * ncv, 0.0 );
    DD_vector z( ldz * ncv, 0.0 );
    DD_vector resid( N, 0.0 );
    DD_vector cfac( lda * N * 2, 0.0 );
    // vector storage of the eigenvalues
    DD_vector dr( nev + 1, 0.0 );
    DD_vector di( nev + 1, 0.0 );
    // number of lower bands ( upper is wired to be the same )
    const unsigned kl( p_A -> noffdiag() );
    // info = 0 => random initial vector
    int info( 0 );
    // max number of iterations
    iparam[ 2 ] = 300;
    // shift and invert mode
    iparam[ 6 ] = 3;
    // the shift value
    double sigmar( shift.real() );
    double sigmai( shift.imag() );
    // convert the banded matrix to LAPACK contiguous FORTRAN format
    FortranData A( *p_A );
    FortranData M( *p_B );
    // call the ARPACK banded wrapper with LM = "largest magnitude"
    ARPACK_DNBAND( &dr[ 0 ], &di[ 0 ], &z[ 0 ], ldz, sigmar,
                   sigmai, N, A.base(), M.base(), lda, &rfac[ 0 ], &cfac[ 0 ], kl,
                   ( char* ) "LM", nev, &resid[ 0 ], ncv, &v[ 0 ], ldv, &iparam[ 0 ],
                   info );
    // number of converged eigenvalues
    const unsigned nconverged( iparam[ 4 ] );
#ifdef DEBUG

    std::cout << "Number of converged eigenvectors = " << nconverged << "\n";
    std::cout << "Requested number of values = " << nev << "\n";
#endif
    // clear out the eigenvalues vector
    all_eigenvalues.clear();
    // store the eigenvalues
    for ( std::size_t i = 0; i < nconverged; ++i )
    {
      all_eigenvalues.push_back( D_complex( dr[ i ], di[ i ] ) );
    }
    if ( calc_eigenvectors )
    {
      /// \todo How are complex conjugate pairs returned
      /// in the ARPACK<double> banded solver?
      /// This is probably broken for complex eigenvalues.

      // clear out the eigenvectors
      all_eigenvectors = CD_matrix( nconverged, N, 0.0 );
      // store the eigenvectors
      for ( std::size_t i = 0; i < nconverged; ++i )
      {
        CD_vector row;
        for ( std::size_t j = 0; j < N; ++j )
        {
          row.push_back( v[ i * N + j ] );
        }
        all_eigenvectors[ i ] = row;
      }
    }
#endif

  }

  template <>
  void BandedLinearEigenSystem< std::complex<double> >::eigensolve_arpack()
  {
#ifndef ARPACK
    std::string problem;
    problem = "The BandedLinearEigenSystem::eigensolve_arpack_without_vectors method has been called\n";
    problem += "but the compiler option -DARPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal( problem );
#else
    // size of the system
    const unsigned N( p_A -> nrows() );
    // number of vectors to be used
    const unsigned ncv( NARNOLDI );
    if ( ( nev + 2 > NARNOLDI ) || ( nev > N ) )
    {
      std::string problem;
      problem = "The BandedLinearEigenSystem::eigensolve_arpack_without_vectors has been called\n";
      problem += "with inconsistent choices for the number of eigenvalues to find and/or \n";
      problem += "the NARNOLDI vectors.\n";
      throw ExceptionExternal( problem, 0 );
    }
    // extra storage factor for complex contents
    const unsigned cmp( 2 );
    // leading dimension of the banded matrix
    const unsigned lda( 3 * p_A -> noffdiag() + 1 );
    // leading dimension of additional work arrays
    const unsigned ldv = N;
    const unsigned ldz = N;
    // parameter vector
    std::vector<int> iparam( 11, 0 );
    // workspaces
    DD_vector fac( cmp * lda * N, 0.0 );
    DD_vector v( cmp * ldv * ncv, 0.0 );
    DD_vector z( cmp * ldz * ncv, 0.0 );
    DD_vector resid( cmp * N, 0.0 );
    // vector storage of the eigenvalues
    DD_vector d( cmp * ( nev + 1 ), 0.0 );
    // number of lower bands ( upper is wired to be the same )
    const unsigned kl( p_A -> noffdiag() );
    // info = 0 => random initial vector
    int info( 0 );
    // max number of iterations
    iparam[ 2 ] = 300;
    // shift and invert mode
    iparam[ 6 ] = 3;
    // the shift value
    DD_vector sigma( cmp, 0.0 );
    sigma[ 0 ] = shift.real();
    sigma[ 1 ] = shift.imag();
    // convert the banded matrix to LAPACK contiguous FORTRAN format
    FortranData A( *p_A );
    FortranData M( *p_B );
    // call the ARPACK banded wrapper with LM = "largest magnitude"
    ARPACK_ZNBAND( &d[ 0 ], &z[ 0 ], ldz, &sigma[ 0 ],
                   N, A.base(), M.base(), lda, &fac[ 0 ], kl,
                   ( char* ) "LM", nev, &resid[ 0 ], ncv, &v[ 0 ], ldv, &iparam[ 0 ],
                   info );
    // number of converged eigenvalues
    const unsigned nconverged( iparam[ 4 ] );
#ifdef DEBUG

    std::cout << "Number of converged eigenvectors = " << nconverged << "\n";
    std::cout << "Requested number of values = " << nev << "\n";
#endif
    // clear out the eigenvalues vector
    all_eigenvalues.clear();
    // store the eigenvalues
    for ( std::size_t i = 0; i < nconverged; i += 2 )
    {
      all_eigenvalues.push_back( D_complex( d[ i ], d[ i + 1 ] ) );
    }
    if ( calc_eigenvectors )
    {
      /// \todo ARPACK<complex> eigenvector retrieval is untested
      /// for the banded solver.

      // clear out the eigenvectors
      all_eigenvectors = CD_matrix( nconverged, N, 0.0 );
      // store the eigenvectors
      for ( std::size_t i = 0; i < nconverged; ++i )
      {
        CD_vector row;
        for ( std::size_t j = 0; j < cmp * N; j += 2 )
        {
          row.push_back( D_complex( v[ i * N * cmp + j ], v[ i * N * cmp + j + 1 ] ) );
        }
        all_eigenvectors[ i ] = row;
      }
    }
#endif

  }


  template class BandedLinearEigenSystem<D_complex>
  ;
  template class BandedLinearEigenSystem<double>
  ;

} // end namespace
