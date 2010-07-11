/// \file DenseLinearEigenSystem.cpp
/// Implementation for the DenseLinearEigenSystem class
/// This class links to LAPACK to perform the solver phase.

#include <vector>
#include <set>

#include <FortranLAPACK.h>
#include <FortranData.h>
#include <DenseLinearEigenSystem.h>
#include <Exceptions.h>
#include <Types.h>

namespace CppNoddy
{

  template <typename _Type>
  DenseLinearEigenSystem<_Type>::DenseLinearEigenSystem( DenseMatrix<_Type>* Aptr, DenseMatrix<_Type>* Bptr ) :
      LinearEigenSystem_base()
  {
    p_A = Aptr;
    p_B = Bptr;
  }

  template <typename _Type>
  DenseLinearEigenSystem<_Type>::~DenseLinearEigenSystem()
  {}

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::eigensolve()
  {
    if ( CALC_EIGENVECTORS )
    {
      eigensolve_lapack_with_vectors();
    }
    else
    {
      eigensolve_lapack_without_vectors();
    }
  }

  template <>
  void DenseLinearEigenSystem<double>::eigensolve_lapack_without_vectors()
  {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal( problem );
#else

    std::size_t N = p_A -> nrows();
    // Cache issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding( 0 );
    if ( ( N % 2 == 0 ) && ( N > 127 ) )
    {
      padding = 1;
    }
#ifdef PARANOID
    if ( ( p_A -> nrows() != p_B -> nrows() ) ||
         ( p_A -> ncols() != p_B -> ncols() ) )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure \n" );
      throw ExceptionGeom( problem, p_A -> nrows(), p_A -> ncols(),
                           p_B -> nrows(), p_B -> ncols() );
    }
#endif
    FortranData Af( *p_A, true, padding );
    FortranData Bf( *p_B, true, padding );
    // eigenvalue storage
    DenseVector<double> alpha_r( N, 0.0 );
    DenseVector<double> alpha_i( N, 0.0 );
    DenseVector<double> beta( N, 0.0 );
    // eigenvector storage
    DenseVector<double> vec_left( 1, 0.0 );
    DenseVector<double> vec_right( 1, 0.0 );
    // some workspace for the LAPACK routine
    DenseVector<double> work( 1, 0.0 );
    int info( 0 );
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_DGGEV( ( char* ) "N", ( char* ) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], -1, info );
    int required_workspace = ( int ) work[ 0 ];
#ifdef DEBUG

    std::cout << " [DEBUG] DenseLinearEigenSystem::eigensolve_lapack_without_vectors is requesting \n";
    std::cout << " [DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize( required_workspace );
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_DGGEV( ( char* ) "N", ( char* ) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], required_workspace, info );
    if ( 0 != info )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure. \n" );
      throw ExceptionExternal( problem, info );
    }
    // create a complex eigenvalue vector
    EIGENVALUES_ALPHA = DenseVector<D_complex>( N, 0.0 );
    for ( std::size_t i = 0; i < N; ++i )
    {
      const D_complex eye( 0.0, 1.0 );
      EIGENVALUES_ALPHA[ i ] = alpha_r[ i ] + alpha_i[ i ] * eye;
    }
    // set the eigenvalue member data
    EIGENVALUES_BETA = beta;
#endif

  }

  template <>
  void DenseLinearEigenSystem< std::complex<double> >::eigensolve_lapack_without_vectors()
  {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal( problem );
#else

    std::size_t N = p_A -> nrows();
    // Cache issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding( 0 );
    if ( ( N % 2 == 0 ) && ( N > 127 ) )
    {
      padding = 1;
    }
#ifdef PARANOID
    if ( ( p_A -> nrows() != p_B -> nrows() ) ||
         ( p_A -> ncols() != p_B -> ncols() ) )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure. \n" );
      throw ExceptionGeom( problem, p_A -> nrows(), p_A -> ncols(),
                           p_B -> nrows(), p_B -> ncols() );
    }
#endif
    // transpose the input matrices so they are in column_major format
    FortranData Af( *p_A, true, padding );
    FortranData Bf( *p_B, true, padding );
    // eigenvalue storage
    DenseVector<double> alpha( 2 * N, 0.0 );
    DenseVector<double> beta( 2 * N, 0.0 );
    // eigenvector storage
    DenseVector<double> vec_left( 2, 0.0 );
    DenseVector<double> vec_right( 2, 0.0 );
    // new the eigenvector storage
    // some workspace for the LAPACK routine
    DenseVector<double> work( 2, 0.0 );
    DenseVector<double> rwork( 8 * N, 0.0 );
    int info( 0 );
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_ZGGEV( ( char* ) "N", ( char* ) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], -1, &rwork[ 0 ], info );
    int required_workspace = int( work[ 0 ] );
#ifdef DEBUG

    std::cout << "[DEBUG] DenseLinearEigenSystem::eigensolve_lapack_without_vectors is requesting \n";
    std::cout << "[DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize( 2 * required_workspace );
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_ZGGEV( ( char* ) "N", ( char* ) "N", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], 1, &work[ 0 ], required_workspace, &rwork[ 0 ], info );
    // error reporting
    if ( 0 != info )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_without_vectors method has detected a failure \n" );
      throw ExceptionExternal( problem, info );
    }
    // create a complex eigenvalue vector for returning
    EIGENVALUES_ALPHA = DenseVector<D_complex>( N, 0.0 );
    EIGENVALUES_BETA = DenseVector<D_complex>( N, 0.0 );
    {
      const D_complex eye( 0.0, 1.0 );
      for ( std::size_t i = 0; i < N; ++i )
      {
        EIGENVALUES_ALPHA[ i ] = alpha[ 2 * i ] + alpha[ 2 * i + 1 ] * eye;
        EIGENVALUES_BETA[ i ] = beta[ 2 * i ] + beta[ 2 * i + 1 ] * eye;
      }
    }
#endif

  }

  template <>
  void DenseLinearEigenSystem< double >::eigensolve_lapack_with_vectors()
  {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal( problem );
#else

    std::size_t N = p_A -> nrows();
    // Cache contention issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding( 0 );
    if ( ( N % 2 == 0 ) && ( N > 127 ) )
    {
      padding = 1;
    }
#ifdef PARANOID
    if ( ( p_A -> nrows() != p_B -> nrows() ) ||
         ( p_A -> ncols() != p_B -> ncols() ) )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure. \n" );
      throw ExceptionGeom( problem, p_A -> nrows(), p_A -> ncols(),
                           p_B -> nrows(), p_B -> ncols() );
    }
#endif
    // Convert to fortran data  incl. a transpose of the
    // input matrices so they are in column_major format then include padding
    FortranData Af( *p_A, true, padding );
    FortranData Bf( *p_B, true, padding );
    // eigenvalue storage
    DenseVector<double> alpha_r( N, 0.0 );
    DenseVector<double> alpha_i( N, 0.0 );
    DenseVector<double> beta( N, 0.0 );
    // new the right eigenvector storage
    DenseVector<double> vec_left( 1, 0.0 );
    DenseVector<double> vec_right( N * N, 0.0 );
    // some workspace for the LAPACK routine
    DenseVector<double> work( 1, 0.0 );
    // return integer for LAPACK
    int info( 0 );
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_DGGEV( ( char* ) "N", ( char* ) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], -1, info );
    int required_workspace = 4 * int( work[ 0 ] );
#ifdef DEBUG

    std::cout << "[DEBUG] DenseLinearEigenSystem::eigensolve_lapack_with_vectors is requesting \n";
    std::cout << "[DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize( required_workspace );
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_DGGEV( ( char* ) "N", ( char* ) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha_r[ 0 ], &alpha_i[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], required_workspace, info );
    if ( 0 != info )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure.\n" );
      throw ExceptionExternal( problem, info );
    }
    // create a complex eigenvalue vector
    EIGENVALUES_ALPHA = DenseVector<D_complex>( N, 0.0 );
    // complex eigenvector matrix
    ALL_EIGENVECTORS = DenseMatrix<D_complex>( N, N, 0.0 );
    // step through the eigenvalues
    for ( std::size_t i = 0; i < N; ++i )
    {
      const D_complex eye( 0.0, 1.0 );
      // make the complex vector of alpha
      EIGENVALUES_ALPHA[ i ] = alpha_r[ i ] + alpha_i[ i ] * eye;
      if ( std::abs( alpha_i[ i ] ) > 0.0 )
      {
        // eigenvector is complex
        for ( std::size_t k = 0; k < N; ++k )
        {
          ALL_EIGENVECTORS[ i ][ k ] = vec_right[ i * N + k ] + eye * vec_right[ ( i + 1 ) * N + k ];
          // store the conjugate too for completeness
          ALL_EIGENVECTORS[ i + 1 ][ k ] = vec_right[ i * N + k ] - eye * vec_right[ ( i + 1 ) * N + k ];
        }
        ++i;
      }
      else // eigenvector is real
      {
        for ( std::size_t k = 0; k < N; ++k )
        {
          ALL_EIGENVECTORS( i, k ) = vec_right[ i * N + k ];
        }
      }
    }
    // set the eigenvalue member data
    EIGENVALUES_BETA = beta;
#endif

  }

  template <>
  void DenseLinearEigenSystem< std::complex<double> >::eigensolve_lapack_with_vectors()
  {
#ifndef LAPACK
    std::string problem;
    problem = "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built.";
    throw ExceptionExternal( problem );
#else

    std::size_t N = p_A -> nrows();
    // Cache issues of varying significance plague problems of size 2^j + 2^k + ...
    // when LDA = N, so this is my shameless 'tweak' to maintain predictable
    // performance, at least for N <=1024 or so.
    int padding( 0 );
    if ( ( N % 2 == 0 ) && ( N > 127 ) )
    {
      padding = 1;
    }
#ifdef PARANOID
    if ( ( p_A -> nrows() != p_B -> nrows() ) ||
         ( p_A -> ncols() != p_B -> ncols() ) )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure. \n" );
      throw ExceptionGeom( problem, p_A -> nrows(), p_A -> ncols(),
                           p_B -> nrows(), p_B -> ncols() );
    }
#endif
    // transpose the input matrices so they are in column_major format
    FortranData Af( *p_A, true, padding );
    FortranData Bf( *p_B, true, padding );
    // eigenvalue storage
    DenseVector<double> alpha( 2 * N, 0.0 );
    DenseVector<double> beta( 2 * N, 0.0 );
    // eigenvector storage
    DenseVector<double> vec_left( 2, 0.0 );
    DenseVector<double> vec_right( 2 * N * N, 0.0 );
    // some workspace for the LAPACK routine
    DenseVector<double> work( 2, 0.0 );
    DenseVector<double> rwork( 8 * N, 0.0 );
    int info( 0 );
    // Call FORTRAN LAPACK to get the required workspace
    LAPACK_ZGGEV( ( char* ) "N", ( char* ) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], -1, &rwork[ 0 ], info );
    int required_workspace = int( work[ 0 ] );
#ifdef DEBUG

    std::cout << "[DEBUG] DenseLinearEigenSystem::eigensolve_lapack_with_vectors is requesting \n";
    std::cout << "[DEBUG] a workspace vector of size " << required_workspace << "\n";
#endif

    work.resize( 2 * required_workspace );
    // call FORTRAN LAPACK again with the optimum workspace
    LAPACK_ZGGEV( ( char* ) "N", ( char* ) "V", N, Af.base(), N + padding, Bf.base(), N + padding, &alpha[ 0 ], &beta[ 0 ], &vec_left[ 0 ], 1, &vec_right[ 0 ], N, &work[ 0 ], required_workspace, &rwork[ 0 ], info );
    if ( 0 != info )
    {
      std::string problem( "The DenseLinearEigenSystem::eigensolve_lapack_with_vectors method has detected a failure.\n" );
      throw ExceptionExternal( problem, info );
    }
    // create a complex eigenvalue vector
    EIGENVALUES_ALPHA = DenseVector<D_complex>( N, 0.0 );
    EIGENVALUES_BETA = DenseVector<D_complex>( N, 0.0 );
    // complex eigenvector matrix
    ALL_EIGENVECTORS = DenseMatrix<D_complex>( N, N, 0.0 );
    // step through the eigenvalues
    for ( std::size_t i = 0; i < N; ++i )
    {
      const D_complex eye( 0.0, 1.0 );
      EIGENVALUES_ALPHA[ i ] = alpha[ 2 * i ] + alpha[ 2 * i + 1 ] * eye;
      EIGENVALUES_BETA[ i ] = beta[ 2 * i ] + beta[ 2 * i + 1 ] * eye;
      for ( std::size_t j = 0; j < N; ++j )
      {
        ALL_EIGENVECTORS( i, j ) = vec_right[ 2 * i * N + 2 * j ] + vec_right[ 2 * i * N + 2 * j + 1 ] * eye;
      }
    }
#endif

  }

  template <typename _Type>
  DenseVector<D_complex> DenseLinearEigenSystem<_Type>::get_tagged_eigenvalues() const
  {
    if ( TAGGED_INDICES.size() == 0 )
    {
      std::string problem;
      problem = "In DenseLinearEigenSystem.get_tagged_eigenvalues() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime( problem );
    }
    // storage for the eigenvalues
    DenseVector<D_complex> evals;
    // loop through the tagged set
    for ( iter p = TAGGED_INDICES.begin(); p != TAGGED_INDICES.end(); ++p )
    {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      // work out the complex eigenvalue associated with this index
      // and add it to the vector
      evals.push_back( EIGENVALUES_ALPHA[ j ] / EIGENVALUES_BETA[ j ] );
    }
    // return the complex vector of eigenvalues
    return evals;
  }

  template <typename _Type>
  DenseMatrix<D_complex> DenseLinearEigenSystem<_Type>::get_tagged_eigenvectors() const
  {
    if ( TAGGED_INDICES.size() == 0 )
    {
      std::string problem;
      problem = "In DenseLinearEigenSystem.get_tagged_eigenvectors() : there are\n";
      problem += "no eigenvalues that have been tagged. This set is empty.\n";
      throw ExceptionRuntime( problem );
    }
    // order of the problem
    std::size_t N = EIGENVALUES_ALPHA.size();
    // eigenvector storage : size() eigenvectors each of length N
    DenseMatrix<D_complex> evecs( TAGGED_INDICES.size(), N, 0.0 );
    std::size_t row = 0;
    // loop through the tagged set
    for ( iter p = TAGGED_INDICES.begin(); p != TAGGED_INDICES.end(); ++p )
    {
      // get the index of the relevant eigenvalue from the set
      std::size_t j = *p;
      // put the eigenvector in the matrix
      evecs[ row ] = ALL_EIGENVECTORS[ j ];
      // next row/eigenvector
      ++row;
    }
    return evecs;
  }

  // EIGENVALUE/VECTOR TAGGING

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_disc( const int &val, const double& radius )
  {
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < EIGENVALUES_ALPHA.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at SHIFT then include it
      if ( std::abs( EIGENVALUES_ALPHA[ i ] - SHIFT * EIGENVALUES_BETA[ i ] ) <
           std::abs( radius * EIGENVALUES_BETA[ i ] ) )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_right( const int &val )
  {
    double real_value = SHIFT.real();
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < EIGENVALUES_ALPHA.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at z then include it
      if ( ( EIGENVALUES_ALPHA[ i ] * std::conj( EIGENVALUES_BETA[ i ] ) ).real() >
           std::pow( std::abs( EIGENVALUES_BETA[ i ] ), 2 ) * real_value )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_left( const int &val )
  {
    double real_value = SHIFT.real();
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < EIGENVALUES_ALPHA.size(); ++i )
    {
      // if the eigenvalue is in the disc centred at z then include it
      if ( ( EIGENVALUES_ALPHA[ i ] * std::conj( EIGENVALUES_BETA[ i ] ) ).real() <
           std::pow( std::abs( EIGENVALUES_BETA[ i ] ), 2 ) * real_value )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_upper( const int &val )
  {
    double imag_value = SHIFT.imag();
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < EIGENVALUES_ALPHA.size(); ++i )
    {
      // if the eigenvalue is in the half plane then include it
      if ( ( EIGENVALUES_ALPHA[ i ] * std::conj( EIGENVALUES_BETA[ i ] ) ).imag() >
           std::pow( std::abs( EIGENVALUES_BETA[ i ] ), 2 ) * imag_value )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }

  template <typename _Type>
  void DenseLinearEigenSystem<_Type>::tag_eigenvalues_lower( const int &val )
  {
    double imag_value = SHIFT.imag();
    // loop through all the eigenvalues
    for ( std::size_t i = 0; i < EIGENVALUES_ALPHA.size(); ++i )
    {
      // if the eigenvalue is in the half plane then include it
      if ( ( EIGENVALUES_ALPHA[ i ] * std::conj( EIGENVALUES_BETA[ i ] ) ).imag() <
           std::pow( std::abs( EIGENVALUES_BETA[ i ] ), 2 ) * imag_value )
      {
        if ( val > 0 )
        {
          // add it to our set of tagged eigenvalues
          TAGGED_INDICES.insert( TAGGED_INDICES.end(), i );
        }
        else
        {
          // remove it from the set if it exists
          TAGGED_INDICES.erase( i );
        }
      }
    }
  }


  template class DenseLinearEigenSystem<D_complex>
  ;
  template class DenseLinearEigenSystem<double>
  ;

} // end namespace
