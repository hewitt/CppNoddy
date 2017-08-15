/// \file SparseMatrix.cpp
/// Implementation of a SPARSE matrix as
/// an STL Vector of SparseVectors, inheriting from Matrix_base.

#include <string>
#include <complex>

#include <SparseMatrix.h>
#include <Matrix_base.h>
#include <Exceptions.h>


namespace CppNoddy
{

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix( const std::size_t& rows, const std::size_t& cols )
      : Matrix_base<_Type>(), NR( rows ), NC( cols )
  {
    MATRIX.reserve( NR );
    SparseVector<_Type> sparse_row( NC );
    for ( std::size_t i = 0; i < NR; ++i )
    {
      MATRIX.push_back( sparse_row );
    }
  }

  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix( const SparseMatrix<_Type>& source, const std::vector<std::size_t>& source_rows ) :
      Matrix_base<_Type>(), NR( source.nrows() ), NC( source.ncols() )
  {
    MATRIX.reserve( NR );
    for ( std::size_t i = 0; i < NR; ++i )
    {
      MATRIX.push_back( source.get_row( source_rows[i] ) );
    }    
  }
  
  template <typename _Type>
  SparseMatrix<_Type>::SparseMatrix( const SparseMatrix<_Type>& source ) :
      Matrix_base<_Type>( source )
  {
    *this = source;
  }

  template <typename _Type>
  inline SparseMatrix<_Type>& SparseMatrix<_Type>::operator=( const SparseMatrix<_Type>& source )
  {
    if ( this == &source )
      return * this;
    MATRIX = source.MATRIX;
    NR = source.NR;
    NC = source.NC;
    return *this;
  }

  template <typename _Type>
  std::size_t SparseMatrix<_Type>::nelts() const
  {
    std::size_t num_of_elts( 0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      num_of_elts += MATRIX[ row ].nelts();
    }
    return num_of_elts;
  }

  template <typename _Type>
  std::size_t SparseMatrix<_Type>::max_in_col( const std::size_t& col,
      const std::size_t& row_min, const std::size_t& row_max ) const
  {
    double maxelt( 0.0 );
    // return outside of the array as default
    std::size_t index = nrows();
    for ( std::size_t row = row_min; row < row_max ; ++row )
    {
      // only bother looking up entries of rows with a first
      // element in a column less than the one we're looking at
      if ( MATRIX[ row ].begin() -> first <= col )
      {
        const double elt( std::abs( MATRIX[ row ].get( col ) ) );
        if ( elt >= maxelt )
        {
          maxelt = elt;
          index = row;
        }
      }
    }
    return index;
  }

  template <typename _Type>
  void SparseMatrix<_Type>::scale( const _Type& mult )
  {
    for ( std::size_t row = 0; row < NR; ++row )
    {
      MATRIX[ row ] *= mult;
    }
  }

  template <typename _Type>
  void SparseMatrix<_Type>::transpose()
  {
    throw ExceptionRuntime( "SparseMatrix.transpose has not been implemented" );
  }

  template <typename _Type>
  double SparseMatrix<_Type>::one_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].one_norm() );
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::two_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].two_norm() );
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::inf_norm() const
  {
    double max( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      max = std::max( max, MATRIX[ row ].inf_norm() );
    }
    return max;
  }

  template <typename _Type>
  double SparseMatrix<_Type>::frob_norm() const
  {
    double sum ( 0.0 );
    for ( std::size_t row = 0; row < NR; ++row )
    {
      sum += MATRIX[ row ].two_norm();
    }
    return sum;
  }

  template <typename _Type>
  DenseVector<_Type> SparseMatrix<_Type>::multiply( const DenseVector<_Type>& X ) const
  {
    throw ExceptionRuntime( "SparseMatrix.multiply has not been implemented" );
  }

#ifdef SUPERLU
  template <typename _Type>
  void SparseMatrix<_Type>::get_row_compressed_superlu( _Type* storage, int* cols, int* rows )
  {
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    //std::size_t last_col( NC + 1 ); // no column should have NC + 1, so this is a dummy start value
    std::size_t i( 0 ); // where we are in the storage vector
    //
    for ( std::size_t row = 0; row < NR; ++row )
    {
      // flag to indicate that we're on a new coloumn
      bool new_row( true );
      pos = MATRIX[ row ].begin();
      do
      {
        _Type elt( pos -> second );
        int col( pos -> first );
        storage[ i ] = elt;
        cols[ i ] = col;
        ++pos;
        if ( new_row )
        {
          rows[ row ] = i;
          new_row = false;
        }
        ++i;
      }
      while ( pos != MATRIX[ row ].end() );
    }
    // last entry points to end + 1
    rows[ NR ] = nelts();
  }
#endif

#ifdef MUMPS_SEQ
  template <>
  void SparseMatrix<double>::get_row_compressed_mumps_seq( double* storage, int* cols, int* rows )
  {
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    std::size_t i( 0 ); // where we are in the storage vector
    //
    for ( std::size_t row = 0; row < NR; ++row )
    {
      // start at the begining of this row
      pos = MATRIX[ row ].begin();
      do
      {
        // for each non-zero elt in the row
        double elt( pos -> second );
        int col( pos -> first );
        storage[ i ] = elt;
        // +1 to return FORTRAN indexing
        cols[ i ] = col+1;
        rows[ i ] = row+1;
        ++pos;
        ++i;
      }
      while ( pos != MATRIX[ row ].end() );
    }
  }


  template <>
  void SparseMatrix<std::complex<double> >::get_row_compressed_mumps_seq( mumps_double_complex* storage, int* cols, int* rows )
  {
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    std::size_t i( 0 ); // where we are in the storage vector
    //
    for ( std::size_t row = 0; row < NR; ++row )
    {
      // start at the begining of this row
      pos = MATRIX[ row ].begin();
      do
      {
        // for each non-zero elt in the row
        mumps_double_complex elt;
        elt.r = real(pos -> second);
        elt.i = imag(pos -> second);
        int col( pos -> first );
        storage[ i ] = elt;
        // +1 to return FORTRAN indexing
        cols[ i ] = col+1;
        rows[ i ] = row+1;
        ++pos;
        ++i;
      }
      while ( pos != MATRIX[ row ].end() );
    }
  }
#endif


#ifdef PETSC_Z
  template <>
  void SparseMatrix<std::complex<double> >::get_row_compressed_petsc( PetscScalar* storage, PetscInt* cols, PetscInt* rows )
  {
    //std::string problem;
    //problem = "The SparseMatrix::get_row_compressed_petsc method was called for a SparseMatrix<D_complex>\n";
    //problem += "even though PETSC_ARCH is currently pointing to a double version of the library.\n";
    //throw ExceptionExternal( problem );
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    std::size_t i( 0 ); // where we are in the storage vector
    //
    for ( std::size_t row = 0; row < NR; ++row )
    {
      // matrix could be singular with an empty row for the mass matrix 
      // of a generalised eigenvalue problem
      if ( MATRIX[row].nelts() > 0 )
      {
        // start at the begining of this row
        pos = MATRIX[ row ].begin();
        do
        {
          // for each non-zero elt in the row
          PetscScalar elt;
          elt = real(pos -> second) + PETSC_i * imag(pos -> second);
          int col( pos -> first );
          storage[ i ] = elt;
          // +1 to return FORTRAN indexing
          cols[ i ] = col;
          rows[ i ] = row;
          ++pos;
          ++i;
        }
        while ( pos != MATRIX[ row ].end() );
      }
    }
  }
#endif
  
#ifdef PETSC_Z
  template <>
  void SparseMatrix<std::complex<double> >::get_row_petsc( PetscInt row, PetscScalar* storage, PetscInt* cols )
  {
    //std::string problem;
    //problem = "The SparseMatrix::get_row_compressed_petsc method was called for a SparseMatrix<D_complex>\n";
    //problem += "even though PETSC_ARCH is currently pointing to a double version of the library.\n";
    //throw ExceptionExternal( problem );
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    std::size_t i(0);
    //
    // matrix could be singular with an empty row for the mass matrix 
    // of a generalised eigenvalue problem
    if ( MATRIX[row].nelts() > 0 )
    {
      // start at the begining of this row
      pos = MATRIX[ row ].begin();
      do
      {
        // for each non-zero elt in the row
        PetscScalar elt;
        elt = std::real(pos -> second) + PETSC_i * std::imag(pos -> second);
        int col( pos -> first );
        storage[ i ] = elt;
        // +1 to return FORTRAN indexing
        cols[ i ] = col;
        ++pos;
        ++i;
      }
      while ( pos != MATRIX[ row ].end() );
    }
  }
#endif
  
#ifdef PETSC_D
  template <>
  void SparseMatrix<double >::get_row_compressed_petsc( PetscScalar* storage, PetscInt* cols, PetscInt* rows )
  {
    //std::string problem;
    //problem = "The SparseMatrix::get_row_compressed_petsc method was called for a SparseMatrix<double>\n";
    //problem += "even though PETSC_ARCH is currently pointing to a complex version of the library.\n";
    //throw ExceptionExternal( problem );
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    std::size_t i( 0 ); // where we are in the storage vector
    //
    for ( std::size_t row = 0; row < NR; ++row )
    {
      // matrix could be singular with an empty row for the mass matrix 
      // of a generalised eigenvalue problem
      if ( MATRIX[row].nelts() > 0 )
      {
        // start at the begining of this row
        pos = MATRIX[ row ].begin();
        do
        {
          // for each non-zero elt in the row
          PetscScalar elt;
          elt = pos -> second;
          int col( pos -> first );
          storage[ i ] = elt;
          // +1 to return FORTRAN indexing
          cols[ i ] = col;
          rows[ i ] = row;
          ++pos;
          ++i;
        }
        while ( pos != MATRIX[ row ].end() );
      }
    }
  }
#endif
  
#ifdef PETSC_D
  template <>
  void SparseMatrix<double >::get_row_petsc( PetscInt row, PetscScalar* storage, PetscInt* cols )
    {
    //std::string problem;
    //problem = "The SparseMatrix::get_row_petsc method was called for a SparseMatrix<double>\n";
    //problem += "even though PETSC_ARCH is currently pointing to a complex version of the library.\n";
    //throw ExceptionExternal( problem );
    // iterator to the maps that are used in SparseVector
    // this is bad form as it exposes the internals of the SparseVector storage
    citer pos;
    std::size_t i(0);
    //
    // matrix could be singular with an empty row for the mass matrix 
    // of a generalised eigenvalue problem
    if ( MATRIX[row].nelts() > 0 )
    {
      // start at the begining of this row
      pos = MATRIX[ row ].begin();
      do
      {
        // for each non-zero elt in the row
        PetscScalar elt;
        elt = pos -> second;
        int col( pos -> first );
        storage[ i ] = elt;
        // +1 to return FORTRAN indexing
        cols[ i ] = col;
        ++pos;
        ++i;
      }
      while ( pos != MATRIX[ row ].end() );
    }
  }
#endif


  template <typename _Type>
  void SparseMatrix<_Type>::dump() const
  {
    std::cout << "SPARSE mtx size = " << NR << " x  sparse \n";
    std::cout.precision( 4 );
    std::cout << "- start matrix \n";
    for ( std::size_t row = 0; row < NR; ++row )
    {
      std::cout << " row " << row << " :  ";
      MATRIX[ row ].dump();
      std::cout << " \n";
    }
    std::cout << "- end matrix \n";
  }

  // the versions to be used are:
  template class SparseMatrix<double>
  ;
  template class SparseMatrix<std::complex<double> >
  ;

} // end namespace

