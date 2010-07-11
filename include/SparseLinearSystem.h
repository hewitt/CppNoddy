/// \file SparseLinearSystem.h
/// Specification of the linear system class.

#ifndef SPARSELINEARSYSTEM_H
#define SPARSELINEARSYSTEM_H

#include <SparseMatrix.h>
#include <DenseVector.h>
#include <Exceptions.h>
#include <LinearSystem_base.h>
#ifdef SUPERLU
#include <slu_ddefs.h>
#endif

namespace CppNoddy
{

  /// A linear system class for vector right-hand sides.
  /// The class is constructed for SPARSE typed problems of the form
  /// \f[ A_{NxN} \,{\underline x}_i = B_{1xN} \f].
  /// Not really useful, just for testing.
  template <typename _Type>
  class SparseLinearSystem : public LinearSystem_base
  {

  public:

    /// Constructor for a sparse linear system object.
    /// \param Aptr A pointer to the 'A matrix', an NxN double/complex sparse matrix
    /// \param Bptr A pointer to the 'B vector' a size N double/complex dense vector
    /// \param which A string that indicates which solver to use
    SparseLinearSystem( SparseMatrix<_Type>* Aptr, DenseVector<_Type>* Bptr, std::string which = "native" );

    /// Destructor for a linear system object.
    ~SparseLinearSystem()
    {}

    /// Solve the sparse system
    void solve();

  private:

    /// Solve the linear system by linking to the SuperLU library
    void solve_superlu()
    {
#ifndef SUPERLU
      std::string problem = "The SparseLinearSystem::solve_superlu method has been called\n";
      problem += "but the compiler option -DSUPERLU was not provided when\n";
      problem += "the library was built and so SuperLU support is not available.";
      throw ExceptionExternal( problem );
#else
      // standard SuperLU preamble straight from the example problems
      SuperMatrix sA, L, U, sB;
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
      p_A -> get_row_compressed( storage, cols, rows );

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
      mem_usage_t   mem_usage;
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

    /// Solve the linear system using the native elimination
    void solve_native();

    /// A wrapped up less_than check that throws an exception
    /// if the matrix elements are complex.
    bool lt( _Type value ) const;

    /// Back substitution routine for dense systems.
    /// \param A The upper triangular matrix LHS
    /// \param B The dense vector RHS
    void backsub( SparseMatrix<_Type> &A, DenseVector<_Type> &B ) const;

    /// check on pivot size for the native elimination routine
    double MIN_PIV;

    /// pointer to a sparse LHS matrix
    SparseMatrix<_Type>* p_A;
    /// pointer to the RHS vector
    DenseVector<_Type>* p_B;
    // boolean flag to indicate if we have the matrix in SuperLU (compressed column) form
    bool have_compressed_array;

  };

} //end namepsace
#endif
