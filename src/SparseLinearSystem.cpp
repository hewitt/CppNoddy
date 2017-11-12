/// \file SparseLinearSystem.cpp
/// Implementation for the LinearSystem class

#include <vector>
#include <set>

#include <SparseLinearSystem.h>
#include <Exceptions.h>
#include <Types.h>
#include <Timer.h>

#ifdef INC_MPI
  #include "mpi.h"
#endif

namespace CppNoddy
{

  template <typename _Type>
  SparseLinearSystem<_Type>::SparseLinearSystem( SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which ) : LinearSystem_base(),
      MIN_PIV( 1.e-12 )
  {
    p_A = Aptr;
    p_B = Bptr;
    VERSION = which;
    //mumps_job_running = false;
    //
    if ( ( VERSION != "superlu" ) && ( VERSION != "native" ) && ( VERSION != "mumps_seq" ) )
    {
      std::string problem;
      problem = "The SparseLinearSystem has been instantiated with an unrecognised\n";
      problem += "request for a solver type. Options: 'native','superlu','mumps_seq'. \n";
      throw ExceptionRuntime( problem );
    }
    #ifdef INC_MPI
      if ( VERSION == "mumps_seq" )
      {
        int flag(0);
        MPI_Initialized( &flag );
        if ( flag != 1 )
        {
          std::string problem;
          problem = "The SparseLinearSystem has been instantiated for a mumps solver.\n";
          problem += "You must run MPI_Init() before calling the mumps solver.\n";
          throw ExceptionRuntime( problem );
        }
      }
    #endif
  }

  template<>
  SparseLinearSystem<double>::~SparseLinearSystem()
  {
    #ifdef MUMPS_SEQ
      if ( mumps_job_running )
      {
        // close things down
        Did_.job = -2;
        dmumps_c(&Did_);
        mumps_job_running = false;
        // mumps_end_job();
        delete[] real_a_;
        delete[] irn_;
        delete[] jcn_;
      }
    #endif
  }

  template<>
  SparseLinearSystem<std::complex<double> >::~SparseLinearSystem()
  {
    #ifdef MUMPS_SEQ
      // std::cout << "[DEBUG] destructor called for a complex SparseLinearSystem with MUMPS_SEQ\n";
      // std::cout << "[DEBUG] mumps_job_running = " << mumps_job_running << "\n";
      if ( mumps_job_running )
      {
        // std::cout << "[DEBUG] mumps_job_running == true\n";
        // close things down
        Zid_.job = -2;
        zmumps_c(&Zid_);
        mumps_job_running = false;
        // mumps_end_job();
        delete[] complex_a_;
        delete[] complex_b_;
        delete[] irn_;
        delete[] jcn_;
      }
    #endif
    // std::cout << "[DEBUG] leaving destructor for complex SparseLinearSystem with MUMPS_SEQ\n";
  }

  template <typename _Type>
  void SparseLinearSystem<_Type>::solve()
  {
    if ( MONITOR_DET )
    {
      std::string problem;
      problem = "You've asked for the LinearSystem object to monitor \n";
      problem += "the determinant for a sparse matrix. This has not been \n";
      problem += "implemented.\n";
      throw ExceptionRuntime( problem );
    }
    #ifdef DEBUG
    Timer timer( "[DEBUG] SparseLinearSystem solver." );
    timer.start();
    #endif
    if ( "superlu" == VERSION )
    {
      solve_superlu();
    }
    else // we catch incorrect VERSION choices in the ctor
    {
      #ifdef MUMPS_SEQ
        if ( "mumps_seq" == VERSION )
        {
          solve_mumps_analysis();
          solve_mumps_factorize();
          solve_mumps_solve_using_factorization();
        }
        else
        {
          solve_native();
        }
      #else
        std::string problem;
        problem = "You've asked for the LinearSystem object to solve \n";
        problem += "using the MUMPS library. This has not been \n";
        problem += "enabled via -DMUMPS_SEQ.\n";
        throw ExceptionRuntime( problem );
      #endif // mumps
    }
    #ifdef DEBUG
    timer.stop();
    timer.print();
    #endif
  }



  template <typename _Type>
  void SparseLinearSystem<_Type>::solve_native()
  {
    const std::size_t Nr = p_A -> nrows();
    // step through rows
    for ( std::size_t l = 0 ; l < Nr - 1 ; ++l )
    {
      // find max index in column 'l' in the range [l, Nr)
      const std::size_t index = p_A -> max_in_col( l , l , Nr );

      // if index of max elt is not the diagonal then swap the rows
      if ( l != index )
      {
        // switch the index row (with maxel) to current position "l"
        p_A -> row_swap( l, index );
        // switch the elts in RHS R-vector
        p_B -> swap( l, index );
      }

      // The diagonal entry is the first :
      // const _Type diag_entry = matrix[ l ].get( l );
      const _Type diag_entry = ( p_A -> MATRIX[ l ] ).begin() -> second;
#ifdef PARANOID
      if ( std::abs( diag_entry ) < MIN_PIV )
      {
        std::string problem( "The pivot in SparseLinearSystem is under the minimum tolerance.\n" );
        throw ExceptionRuntime( problem );
      }
#endif
      // eliminate all entries below
      for ( std::size_t row = l + 1 ; row < Nr ; ++row )
      {
        // subtract rows: R_row = R_row - mult * R_l
        // but optimise out the zero elements to the left
        // cycle through R_l data
        typename std::map< std::size_t, _Type >::const_iterator pos_ro = ( p_A -> MATRIX[ l ] ).begin();
        typename std::map< std::size_t, _Type >::iterator pos_rw = ( p_A -> MATRIX[ row ] ).begin();
        if ( pos_rw -> first <= l )
        {
          // this should be equivalent to 'MATRIX[ row ].get( l ) / diag_entry;"
          const _Type mult = ( pos_rw -> second ) / diag_entry;
          do
          {
            std::size_t index_rw = pos_rw -> first;
            std::size_t index_ro = pos_ro -> first;
            if ( index_rw == index_ro )
            {
              // element in both vectors
              if ( pos_rw -> first == l )
              {
                // this entry should now be eliminated from the matrix
                ( p_A -> MATRIX[ row ] ).erase( pos_rw );
              }
              else
              {
                pos_rw -> second -= ( pos_ro -> second ) * mult;
              }
              ++pos_rw;
              ++pos_ro;
            }
            if ( index_rw > index_ro )
            {
              // element is in X but not 'this'
              ( p_A -> MATRIX[ row ] ).set( index_ro ) = -( pos_ro -> second ) * mult;
              ++pos_ro;
            }
            if ( index_rw < index_ro )
            {
              // element is in 'this' but not X
              ++pos_rw;
            }
          }
          while ( pos_ro != ( p_A -> MATRIX[ l ] ).end() &&
                  pos_rw != ( p_A -> MATRIX[ row ] ).end() );

          if ( pos_ro != ( p_A -> MATRIX[ l ] ).end() )
          {
            // need to finish the X data
            do
            {
              ( p_A -> MATRIX[ row ] ).set( pos_ro -> first ) = -( pos_ro -> second ) * mult;
              ++pos_ro;
            }
            while ( pos_ro != ( p_A -> MATRIX[ l ] ).end() );
          }
          // this is a scalar operation
          p_B -> operator[] ( row ) -= p_B -> operator[] ( l ) * mult;
        }
      }

    } // close l-loop
#ifdef PARANOID
    // check last row for singular matrix
    // The diagonal entry should be the first:
    // const _Type diag_entry = matrix[ Nr - 1 ].get( Nr - 1 );
    const _Type diag_entry = ( p_A -> MATRIX[ Nr - 1 ] ).begin() -> second;
    if ( std::abs( diag_entry ) < MIN_PIV )
    {
      std::string problem( "The pivot in NSMatrix.GJE is under the minimum tolerance.\n" );
      throw ExceptionRuntime( problem );
    }
#endif
    backsub( *p_A, *p_B );
  }




  template <>
  void SparseLinearSystem<double>::solve_superlu()
  {
#ifndef SUPERLU
    std::string problem = "The SparseLinearSystem::solve_superlu method has been called\n";
    problem += "but the compiler option -DSUPERLU was not provided when\n";
    problem += "the library was built and so SuperLU support is not available.";
    throw ExceptionExternal( problem );
#else
    using namespace SLUD;
    // standard SuperLU preamble straight from the example problems
    SuperMatrix sA, L, U, sB;
    // All these integers are 4 bytes (32bit) -- presumably because the
    // superlu library is compiled in this way, as is BLAS probably.
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int info;
    int m = p_A -> nrows();
    int n = p_A -> ncols();
    int nnz = p_A -> nelts();
    int nrhs = 1;
    superlu_options_t options;
    SuperLUStat_t stat;

    // convert our SparseMatrix problem into compressed column data
    double* storage;
    int* rows;
    int* cols;
    if ( !( storage = doubleMalloc( nnz ) ) )
      ABORT( "Malloc fails for storage[]." );
    if ( !( cols = intMalloc( nnz ) ) )
      ABORT( "Malloc fails for rows[]." );
    if ( !( rows = intMalloc( m + 1 ) ) )
      ABORT( "Malloc fails for cols[]." );

    // this is the only intersection with the CppNoddy container
    // this returns all the required row_compressed data
    p_A -> get_row_compressed_superlu( storage, cols, rows );

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix( &sA, m, n, nnz, storage, cols,
                            rows, SLU_NR, SLU_D, SLU_GE );
    // ^ the SLU_NR here indicates that it's row-compressed
    //   the SLU_D indicates double precision
    //   the SLU_GE indicates that its a general matrix with no special properties
    dCreate_Dense_Matrix( &sB, m, nrhs, &( ( *p_B )[0] ), m, SLU_DN, SLU_D, SLU_GE );
    // ^ the SLU_DN indicates that its double & in Fortran column-first format
    //     not that it makes any difference if nrhs=1

    if ( !( perm_r = intMalloc( m ) ) )
      ABORT( "Malloc fails for perm_r[]." );
    if ( !( perm_c = intMalloc( n ) ) )
      ABORT( "Malloc fails for perm_c[]." );

    /* Set the default input options. */
    //options.ColPerm = NATURAL;
    set_default_options( &options );

    /* Initialize the statistics variables. */
    StatInit( &stat );

    // solve & get the solution
    dgssv( &options, &sA, perm_c, perm_r, &L, &U, &sB, &stat, &info );
    double *sol = ( double* ) ( ( DNformat* ) sB.Store ) -> nzval;
    for ( int j = 0; j < n; ++j )
    {
      // return via the CppNoddy DenseVector container
      ( *p_B )[ j ] = sol[ j ];
    }

#ifdef DEBUG
    SCformat *Lstore;
    NCformat *Ustore;
    SLUD::mem_usage_t   mem_usage;
    Lstore = ( SCformat * ) L.Store;
    Ustore = ( NCformat * ) U.Store;
    printf( "[DEBUG] No of nonzeros in factor L = %d\n", Lstore->nnz );
    printf( "[DEBUG] No of nonzeros in factor U = %d\n", Ustore->nnz );
    printf( "[DEBUG] No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n );
    printf( "[DEBUG] FILL ratio = %.1f\n", ( float )( Lstore->nnz + Ustore->nnz - n ) / nnz );
    dQuerySpace( &L, &U, &mem_usage );
    printf( "[DEBUG] L\\U MB %.3f\ttotal MB needed %.3f\n",
            mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6 );
#endif

    SUPERLU_FREE ( perm_r );
    SUPERLU_FREE ( perm_c );
    // usual SuperLU free memory allocated to storage
    Destroy_CompCol_Matrix( &sA );
    Destroy_SuperMatrix_Store( &sB );
    Destroy_SuperNode_Matrix( &L );
    Destroy_CompCol_Matrix( &U );
    StatFree( &stat );
#endif
  }




  template <>
  void SparseLinearSystem<std::complex<double> >::solve_superlu()
  {
#ifndef SUPERLU
    std::string problem = "The SparseLinearSystem::solve_superlu method has been called\n";
    problem += "but the compiler option -DSUPERLU was not provided when\n";
    problem += "the library was built and so SuperLU support is not available.";
    throw ExceptionExternal( problem );
#else
    using namespace SLUZ;
    // standard SuperLU preamble straight from the example problems
    SuperMatrix sA, L, U, sB;
    // All these integers are 4 bytes (32bit) -- presumably because the
    // superlu library is compiled in this way, as is BLAS probably.
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int info;
    int m = p_A -> nrows();
    int n = p_A -> ncols();
    int nnz = p_A -> nelts();
    int nrhs = 1;
    superlu_options_t options;
    SuperLUStat_t stat;

    // convert our SparseMatrix problem into compressed column data
    doublecomplex* storage;
    int* rows;
    int* cols;
    if ( !( storage = doublecomplexMalloc( nnz ) ) )
      ABORT( "Malloc fails for storage[]." );
    if ( !( cols = intMalloc( nnz ) ) )
      ABORT( "Malloc fails for rows[]." );
    if ( !( rows = intMalloc( m + 1 ) ) )
      ABORT( "Malloc fails for cols[]." );

    // temporary row_compressed storage
    //std::complex<double> temp_storage[ nnz ];
    std::vector<std::complex<double> > temp_storage( nnz, 0.0 );
    // this is the only intersection with the CppNoddy container
    // this returns all the required row_compressed data
    p_A -> get_row_compressed_superlu( &(temp_storage[0]), cols, rows );
    // SUPERLU uses a "doublecomplex" struct, so we have to convert
    // to that from the std::complex<double> class for both the
    // matrix and the RHS
    for ( int  k = 0; k < nnz; ++k )
    {
      storage[ k ].r = std::real(temp_storage[ k ]);
      storage[ k ].i = std::imag(temp_storage[ k ]);
    }
    doublecomplex B[ m ];
    for ( int k = 0; k < m; ++k )
    {
      B[ k ].r = real((*p_B)[ k ]);
      B[ k ].i = imag((*p_B)[ k ]);
    }

    /* Create matrix A in the format expected by SuperLU. */
    zCreate_CompCol_Matrix( &sA, m, n, nnz, storage, cols,
                            rows, SLU_NR, SLU_Z, SLU_GE );
    // ^ the SLU_NR here indicates that it's row-compressed
    //   the SLU_D indicates double precision
    //   the SLU_GE indicates that its a general matrix with no special properties
    zCreate_Dense_Matrix( &sB, m, nrhs, &B[0], m, SLU_DN, SLU_Z, SLU_GE );
    // ^ the SLU_DN indicates that its double & in Fortran column-first format
    //     not that it makes any difference if nrhs=1

    if ( !( perm_r = intMalloc( m ) ) )
      ABORT( "Malloc fails for perm_r[]." );
    if ( !( perm_c = intMalloc( n ) ) )
      ABORT( "Malloc fails for perm_c[]." );

    /* Set the default input options. */
    //options.ColPerm = NATURAL;
    set_default_options( &options );

    /* Initialize the statistics variables. */
    StatInit( &stat );

    // solve & get the solution
    zgssv( &options, &sA, perm_c, perm_r, &L, &U, &sB, &stat, &info );
    doublecomplex *sol = ( doublecomplex* ) ( ( DNformat* ) sB.Store ) -> nzval;
    const std::complex<double> eye( 0., 1. );
    for ( int j = 0; j < n; ++j )
    {
      // return via the CppNoddy DenseVector container
      ( *p_B )[ j ] = sol[ j ].r + eye * sol[ j ].i;
    }

#ifdef DEBUG
    SCformat *Lstore;
    NCformat *Ustore;
    SLUZ::mem_usage_t   mem_usage;
    Lstore = ( SCformat * ) L.Store;
    Ustore = ( NCformat * ) U.Store;
    printf( "[DEBUG] No of nonzeros in factor L = %d\n", Lstore->nnz );
    printf( "[DEBUG] No of nonzeros in factor U = %d\n", Ustore->nnz );
    printf( "[DEBUG] No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n );
    printf( "[DEBUG] FILL ratio = %.1f\n", ( float )( Lstore->nnz + Ustore->nnz - n ) / nnz );
    zQuerySpace( &L, &U, &mem_usage );
    printf( "[DEBUG] L\\U MB %.3f\ttotal MB needed %.3f\n",
            mem_usage.for_lu / 1e6, mem_usage.total_needed / 1e6 );
#endif

    SUPERLU_FREE ( perm_r );
    SUPERLU_FREE ( perm_c );
    // usual SuperLU free memory allocated to storage
    Destroy_CompCol_Matrix( &sA );
    Destroy_SuperMatrix_Store( &sB );
    Destroy_SuperNode_Matrix( &L );
    Destroy_CompCol_Matrix( &U );
    StatFree( &stat );
#endif
  }



#ifdef MUMPS_SEQ

//   template <>
//   void SparseLinearSystem<double>::solve_mumps_seq()
//   {
// #ifndef MUMPS_SEQ
//     std::string problem = "The SparseLinearSystem<double>::solve_mumps_seq method has been called\n";
//     problem += "but the compiler option -DMUMPS_SEQ was not provided when\n";
//     problem += "the library was built and so MUMPS_SEQ support is not available.";
//     throw ExceptionExternal( problem );
// #else
//
//     Timer timer;
//
//     DMUMPS_STRUC_C id;
//
//     int ierr, myid;
//     // create a new MPI singleton instance only if none already exists
//     MPIinit* p_library = MPIinit::getInstance();
//     ierr = MPI_Comm_rank(p_library->get_Comm(), &myid);
//     //ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//
//     int nr = p_A -> nrows();
//     int nnz = p_A -> nelts();
//
//     // initialisation is done with -1 as the job
//     id.job=-1;
//     // par has to be set to 1 for sequential or for host to contribute to parallel
//     id.par=1;
//     // not symmetric
//     id.sym=0;
//     // something to handle fortran<->C as specified in MUMPS manual
//     id.comm_fortran=-987654;
//     timer.start();
//     // set things up
//     dmumps_c(&id);
//     timer.stop();
//     std::cout << "MUMPS set up:";
//     timer.print();
//     timer.reset();
//
//     // All these integers are 4 bytes (32bit) -- presumably because the
//     // MUMPS library is compiled in this way, as is BLAS probably.
//     //
//     // we can't do "double a[nnz];" etc b/c this would allocate from the stack and fail
//     // when nnz is not particularly large
//     //
//     // storage for the row indices
//     int* irn = new int[nnz];
//     // storage for the column indices
//     int* jcn = new int[nnz];
//     // storage for the element list
//     double* a = new double[nnz];
//
//     // we don't need alternative storage for the RHS when the problem is real
//     // double rhs[nr];
//     // convert the sparse matrix into the coordinate format required by MUMPS
//     timer.start();
//     p_A -> get_row_compressed_mumps_seq( a, jcn, irn );
//     timer.stop();
//     std::cout << "MUMPS get_row_compressed_mumps_seq():";
//     timer.print();
//     timer.reset();
//
//     if (myid == 0)  // define the problem on the host
//     {
//       // set up the details of the matrix
//       id.n = nr; id.nz =nnz; id.irn=irn; id.jcn=jcn;
//       // set the LHS
//       id.a = a;
//       // keep the RHS in its DenseVector format and just pass the base address to the stored vector
//       id.rhs = &( ( *p_B )[0] );
//     }
//
//     //#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
//     ///* No outputs */
//     //id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
//
//     // these indices are -1 from that in the MUMPS (Fortran indexed) documentation
//     id.icntl[0] = 6; //error output stream (Fortran numbering)
// #ifdef DEBUG
//     id.icntl[1] = 6; //diagnostic/warning/statistics collected output stream (Fortran numbering)
// #else
//     id.icntl[1] = 0; //diagnostic/warning/statistics collected output stream (Fortran numbering)
// #endif
//     id.icntl[2] = 0; //global information collected on host
//     id.icntl[3] = 1; //level of printing verbosity for above streams
//
//     ////ordering metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)
//     //id.icntl[6] = 2; // 7 = default which indicates automatic
//     timer.start();
//     // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
//     id.job = 1;
//     dmumps_c(&id);
//     timer.stop();
//     std::cout << "MUMPS analysis job = 1:";
//     timer.print();
//     timer.reset();
//
//     timer.start();
//     // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
//     id.job = 2;
//     dmumps_c(&id);
//     timer.stop();
//     std::cout << "MUMPS factorization job = 2:";
//     timer.print();
//     timer.reset();
//
//     timer.start();
//     // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
//     id.job = 3;
//     dmumps_c(&id);
//     timer.stop();
//     std::cout << "MUMPS factorization job = 3:";
//     timer.print();
//     timer.reset();
//
//     // if we are interested in timings, we may as well throw in the memory used too
// #ifdef DEBUG
//     std::cout << "[DEBUG] Memory used by MUMPS_SEQ = " << id.info[ 15 ] << "MB.\n";
// #endif
//     // close things down
//     timer.start();
//     // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
//     id.job = -2;
//     dmumps_c(&id);
//     timer.stop();
//     std::cout << "MUMPS finalization job = -2:";
//     timer.print();
//     timer.reset();
//
//     delete[] a;
//     delete[] irn;
//     delete[] jcn;
// #ifdef DEBUG
//     if ( myid == 0 )
//     {
//       std::cout << " [DEBUG] MPI ierr in MUMPS solve phase = " << ierr << "\n";
//     }
// #endif // DEBUG
// #endif // MUMPS ndef
//   }
//
//
//
//   template <>
//   void SparseLinearSystem<D_complex>::solve_mumps_seq()
//   {
// #ifndef MUMPS_SEQ
//     std::string problem = "The SparseLinearSystem<D_complex>::solve_mumps_seq method has been called\n";
//     problem += "but the compiler option -DMUMPS_SEQ was not provided when\n";
//     problem += "the library was built and so MUMPS_SEQ support is not available.";
//     throw ExceptionExternal( problem );
// #else
//     ZMUMPS_STRUC_C id;
//
//     // create a new MPI singleton instance only if none already exists
//     MPIinit* p_library = MPIinit::getInstance();
//     int ierr, myid;
//     ierr = MPI_Comm_rank(p_library->get_Comm(), &myid);
//
//     // All these integers are 4 bytes (32bit) -- presumably because the
//     // MUMPS library is compiled in this way, as is BLAS probably.
//     int nr = p_A -> nrows();
//     int nnz = p_A -> nelts();
//
//     // instance is initialised via job=-1
//     id.job=-1;
//     // for the sequential non-MPI version we need par=1
//     id.par=1;
//     // assumes no symmetry
//     id.sym=0;
//     // something to handle fortran<->C as specified in MUMPS manual
//     id.comm_fortran=-987654;
//     // set things up
//     zmumps_c(&id);
//
//     // storage for row indices
//     int* irn = new int[nnz];
//     // storage for the column indices
//     int* jcn = new int[nnz];
//     // storage for the element list
//     mumps_double_complex* a = new mumps_double_complex[nnz];
//     // RHS vector storage is needed because it is a complex problem
//     mumps_double_complex* rhs = new mumps_double_complex[nr];
//
//     // convert the sparse matrix into the coordinate format required by MUMPS
//     p_A -> get_row_compressed_mumps_seq( a, jcn, irn );
//
//     // set up the RHS vector from the CppNoddy RHS
//     for ( int k = 0; k<nr; ++k )
//     {
//       (rhs[k]).r = real( (*p_B)[k] );
//       (rhs[k]).i = imag( (*p_B)[k] );
//     }
//
//     if (myid == 0)  // define the problem on the host
//     {
//       // set up the details of the matrix
//       id.n = nr; id.nz =nnz; id.irn=irn; id.jcn=jcn;
//       // set the LHS
//       id.a = a;
//       // define the RHS
//       id.rhs = rhs;
//     }
//
//     // these indices are -1 from that in the MUMPS (Fortran indexed) documentation
// #ifdef DEBUG
//     id.icntl[1] = 6; //diagnostic/warning/statistics collected output stream (Fortran numbering)
// #else
//     id.icntl[1] = 0; //diagnostic/warning/statistics collected output stream (Fortran numbering)
// #endif
//     id.icntl[2] = 0; //global information collected on host
//     id.icntl[3] = 1; //level of printing verbosity for above streams
//
//     //ordering metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)
//     id.icntl[6] = 2;
//
//     // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
//     id.job = 6;
//     zmumps_c(&id);
//
//     // if we are interested in timings, we may as well throw in the memory used too
// #ifdef DEBUG
//     std::cout << "[DEBUG] Memory used by MUMPS_SEQ = " << id.info[ 15 ] << "MB.\n";
//     if ( myid == 0 )
//     {
//       std::cout << " [DEBUG] MPI ierr in MUMPS solve phase = " << ierr << "\n";
//     }
// #endif
//
//     const std::complex<double> eye( 0., 1. );
//     for ( std::size_t j = 0; j < nr; ++j )
//     {
//       // return via the CppNoddy DenseVector container
//       ( *p_B )[ j ] = (rhs[ j ]).r + eye * (rhs[ j ]).i;
//     }
//
//     // close down this instance
//     id.job = -2;
//     zmumps_c(&id);
//
//     delete[] a;
//     delete[] rhs;
//     delete[] irn;
//     delete[] jcn;
// #endif
//   }


  template <>
  void SparseLinearSystem<double>::mumps_end_job()
  {
    std::cout << "Closing down real mumps\n";
    // close things down
    Did_.job = -2;
    dmumps_c(&Did_);
    mumps_job_running = false;
    std::cout << "Closed down real mumps\n";
  }




  template <>
  void SparseLinearSystem<std::complex<double> >::mumps_end_job()
  {
    std::cout << "Closing down Complex mumps\n";
    std::cout << "[DEBUG] ending job: mumps_job_running = " << mumps_job_running << "\n";
    // close things down
    Zid_.job = -2;
    zmumps_c(&Zid_);
    mumps_job_running = false;
    std::cout << "Closed down complex mumps\n";
    std::cout << "[DEBUG] ending job: mumps_job_running = " << mumps_job_running << "\n";
  }




  template <>
  void SparseLinearSystem<double>::solve_mumps_analysis()
  {
    // if we have already analyzed/factorized then we need to end
    // that instance and start again to avoid memory leaks
    if ( mumps_job_running )
    {
      mumps_end_job();
    }

    mumps_job_running = true;
    // Initialize the private member data to the appropriate struct
    Did_ = DMUMPS_STRUC_C();

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // analysis must be preceded by initialisation
    // initialisation is done with -1 as the job
    Did_.job=-1;
    // par has to be set to 1 for sequential OR for host to contribute to parallel solve
    Did_.par=1;
    // not symmetric
    Did_.sym=0;
    // something to handle fortran<->C as specified in MUMPS manual
    Did_.comm_fortran=-987654;
    //Did_.comm_fortran = (MUMPS_INT) (p_library->get_Comm());
    // set things up
    dmumps_c(&Did_);

    // All these integers are 4 bytes (32bit) -- presumably because the
    // MUMPS library is compiled in this way, as is BLAS probably.
    //
    //
    int nr = p_A -> nrows();
    int nnz = p_A -> nelts();
    // storage for the row indices
    // we can't do "double a[nnz];" etc b/c this would allocate from the stack and fail
    // when nnz is not particularly large
    irn_ = new int[nnz];
    // storage for the column indices
    jcn_ = new int[nnz];
    // storage for the element list
    real_a_ = new double[nnz];
    // we don't need alternative storage for the RHS when the problem is real
    // as we can just pass the base pointer to a DenseVector
    //
    // convert the sparse matrix into the coordinate format required by MUMPS
    // this is typically a few hundred miliseconds -- so not a big deal
    //p_A -> get_row_compressed_mumps_seq( a_, &jcn_[0], &irn_[0] );
    p_A -> get_row_compressed_mumps_seq( real_a_, jcn_, irn_ );

    if (myid == 0)  // problem must be defined on the host only
    {
      // set up the details of the matrix
      Did_.n = nr; Did_.nz =nnz; Did_.irn=irn_; Did_.jcn=jcn_;
      // set the LHS
      Did_.a = real_a_;
      // keep the RHS in its DenseVector format and just pass the base address to the stored vector
      Did_.rhs = &( ( *p_B )[0] );
    }

    // these indices are -1 from that in the MUMPS (Fortran indexed) documentation
    Did_.icntl[0] = 6; //error output stream (Fortran numbering)
#ifdef DEBUG
    Did_.icntl[1] = 6; //diagnostic/warning/statistics collected output stream (Fortran numbering)
#else
    Did_.icntl[1] = 0; //diagnostic/warning/statistics collected output stream (Fortran numbering)
#endif
    Did_.icntl[2] = 0; //global information collected on host
    Did_.icntl[3] = 1; //level of printing verbosity for above streams
    ////ordering metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)
    Did_.icntl[6] = 5; // 7 = default which indicates automatic

    // do the analysis
    Timer timer("MUMPS analysis phase");
    timer.start();
    Did_.job = 1;
    dmumps_c(&Did_);
    timer.stop();
    timer.print();

    // if we are interested in timings, we may as well throw in the memory used too
#ifdef DEBUG
    std::cout << "[DEBUG] Memory used by MUMPS_SEQ = " << Did_.info[ 15 ] << "MB.\n";
#endif
  }




  template <>
  void SparseLinearSystem<double>::solve_mumps_factorize()
  {
    Timer timer("MUMPS factorization phase");
    timer.start();
    Did_.job = 2;
    dmumps_c(&Did_);
    timer.stop();
    timer.print();
  }





  template <>
  void SparseLinearSystem<double>::solve_mumps_solve_using_factorization()
  {
    // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
    Timer timer("MUMPS solve phase");
    timer.start();
    Did_.job = 3;
    dmumps_c(&Did_);
    timer.stop();
    timer.print();
    #ifdef DEBUG
      std::cout << "[DEBUG] Memory used by MUMPS_SEQ = " << Did_.info[ 15 ] << "MB.\n";
    #endif
  }




  template <>
  void SparseLinearSystem<std::complex<double> >::solve_mumps_analysis()
  {
    // if we have already analyzed/factorized then we need to end
    // that instance and start again to avoid memory leaks
    if ( mumps_job_running )
    {
      mumps_end_job();
    }
    mumps_job_running = true;
    // Initialize the private member data to the appropriate struct
    Zid_ = ZMUMPS_STRUC_C();

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    // analysis must be preceded by initialisation
    // initialisation is done with -1 as the job
    Zid_.job=-1;
    // par has to be set to 1 for sequential OR for host to contribute to parallel solve
    Zid_.par=1;
    // not symmetric
    Zid_.sym=0;
    // something to handle fortran<->C as specified in MUMPS manual
    Zid_.comm_fortran=-987654;
    // set things up
    zmumps_c(&Zid_);

    // All these integers are 4 bytes (32bit) -- presumably because the
    // MUMPS library is compiled in this way, as is BLAS probably.
    //
    //
    int nr = p_A -> nrows();
    int nnz = p_A -> nelts();
    // storage for the row indices
    // we can't do "double a[nnz];" etc b/c this would allocate from the stack and fail
    // when nnz is not particularly large
    irn_ = new int[nnz];
    // storage for the column indices
    jcn_ = new int[nnz];
    // storage for the element list
    complex_a_ = new mumps_double_complex[nnz];
    complex_b_ = new mumps_double_complex[nr];
    // set up the RHS vector from the CppNoddy RHS
    for ( int k = 0; k<nr; ++k )
    {
      (complex_b_[k]).r = real( (*p_B)[k] );
      (complex_b_[k]).i = imag( (*p_B)[k] );
    }

    //

    // convert the sparse matrix into the coordinate format required by MUMPS
    // this is typically a few hundred miliseconds -- so not a big deal
    //p_A -> get_row_compressed_mumps_seq( a_, &jcn_[0], &irn_[0] );
    p_A -> get_row_compressed_mumps_seq( &complex_a_[0], &jcn_[0], &irn_[0] );

    if (myid == 0)  // problem must be defined on the host only
    {
      // set up the details of the matrix
      Zid_.n = nr; Zid_.nz =nnz; Zid_.irn=irn_; Zid_.jcn=jcn_;
      // set the LHS
      Zid_.a = complex_a_;
      // keep the RHS in its DenseVector format and just pass the base address to the stored vector
      Zid_.rhs = complex_b_;
    }

    // these indices are -1 from that in the MUMPS (Fortran indexed) documentation
    Zid_.icntl[0] = 6; //error output stream (Fortran numbering)
#ifdef DEBUG
    Zid_.icntl[1] = 6; //diagnostic/warning/statistics collected output stream (Fortran numbering)
#else
    Zid_.icntl[1] = 0; //diagnostic/warning/statistics collected output stream (Fortran numbering)
#endif
    Zid_.icntl[2] = 0; //global information collected on host
    Zid_.icntl[3] = 1; //level of printing verbosity for above streams
    ////ordering metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)
    Zid_.icntl[6] = 5; // 7 = default which indicates automatic

    // do the analysis
    Timer timer("MUMPS analysis phase");
    timer.start();
    Zid_.job = 1;
    zmumps_c(&Zid_);
    timer.stop();
    timer.print();

    // if we are interested in timings, we may as well throw in the memory used too
#ifdef DEBUG
    std::cout << "[DEBUG] Memory used by MUMPS_SEQ = " << Zid_.info[ 15 ] << "MB.\n";
#endif
  }




  template <>
  void SparseLinearSystem<std::complex<double> >::solve_mumps_factorize()
  {
    Timer timer("MUMPS factorization phase");
    timer.start();
    Zid_.job = 2;
    zmumps_c(&Zid_);
    timer.stop();
    timer.print();
  }





  template <>
  void SparseLinearSystem<std::complex<double> >::solve_mumps_solve_using_factorization()
  {
    // 1 = analysis, 2 = factorization, 3 = solve previously factorised, 4 = same as 1+2, 5 = same as 2+3, 6 = same as 1+2+3
    Timer timer("MUMPS solve phase");
    timer.start();
    Zid_.job = 3;
    zmumps_c(&Zid_);
    const std::complex<double> eye( 0., 1. );
    for ( std::size_t j = 0; j < Zid_.n; ++j )
    {
      // return via the CppNoddy DenseVector container
      ( *p_B )[ j ] = (complex_b_[ j ]).r + eye * (complex_b_[ j ]).i;
    }
    timer.stop();
    timer.print();
    #ifdef DEBUG
      std::cout << "[DEBUG] Memory used by MUMPS_SEQ = " << Zid_.info[ 15 ] << "MB.\n";
    #endif
  }

#endif // mumps


  template <typename _Type>
  void SparseLinearSystem<_Type >::backsub( SparseMatrix<_Type> &A, DenseVector<_Type> &B ) const
  {
    const std::size_t Nr( B.size() );
    DenseVector<_Type> x( Nr, 0.0 );
    // This line should be equivalent to:
    // x[ Nr - 1 ] = B[ Nr - 1 ] / matrix[ Nr - 1 ].get( Nr - 1 );
    x[ Nr - 1 ] = B[ Nr - 1 ] / ( A.MATRIX[ Nr - 1 ].begin() -> second );
    // Note the unusual (aka hacked up) row < N termination condition.
    // We can't do "row <= 0" with std::size_t.
    for ( std::size_t row = Nr - 2; row < Nr; --row )
    {
      _Type sum = 0.0;
      for ( typename std::map< std::size_t, _Type >::const_iterator pos_ro =
              A.MATRIX[ row ].begin(); pos_ro != A.MATRIX[ row ].end(); ++pos_ro )
      {
        std::size_t index_ro = pos_ro -> first;
        if ( index_ro > row )
        {
          sum += ( pos_ro -> second ) * x[ index_ro ];
        }
      }
      // The diagonal entry is the first since it's triangular
      x[ row ] = ( B[ row ] - sum ) / ( A.MATRIX[ row ].begin() -> second );
    }
    B = x;
  }






  template class SparseLinearSystem<D_complex>
  ;
  template class SparseLinearSystem<double>
  ;

} // end namespace
