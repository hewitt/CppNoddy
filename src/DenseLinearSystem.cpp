/// \file DenseLinearSystem.cpp
/// Implementation for the LinearSystem class

#include <vector>

#include <FortranLAPACK.h>
#include <FortranData.h>
#include <DenseLinearSystem.h>
#include <Exceptions.h>

namespace CppNoddy
{

  template <typename _Type>
  DenseLinearSystem<_Type>::DenseLinearSystem( DenseMatrix<_Type>* p_A,
      DenseVector<_Type>* p_B,
      std::string which ) :
      LinearSystem_base(),
      MIN_PIV( 1.e-12 )
  {
    this -> p_A = p_A;
    this -> p_B = p_B;
    VERSION = which;
    if ( ( VERSION != "lapack" ) && ( VERSION != "native" ) )
    {
      std::string problem;
      problem = "The DenseLinearSystem has been instantiated with an unrecognised\n";
      problem += "request for a solver type. Options are 'native' or 'lapack'. \n";
      throw ExceptionRuntime( problem );
    }
  }

  template <typename _Type>
  void DenseLinearSystem<_Type>::solve()
  {
    if ( "lapack" == VERSION )
    {
      solve_lapack();
    }
    else // we catch incorrect VERSION choices in the ctor
    {
      solve_native();
    }
  }

  template <typename _Type>
  void DenseLinearSystem<_Type>::solve_lapack( )
  {
    std::string problem;
    problem = "The solve method for a DenseLinearSystem has not been implemented\n";
    problem += "for the element type used here. \n";
    throw ExceptionExternal( problem );
  }

  // lapack LU solver for double-element dense matrices
  template <>
  void DenseLinearSystem<double>::solve_lapack( )
  {
#ifndef LAPACK
    std::string problem = "The DenseLinearSystem::solve_lapack method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built and so LAPACK support is not available.";
    throw ExceptionExternal( problem );
#else

    std::size_t N = p_A -> nrows();
    FortranData Af( *p_A );
    // a vector to store the pivots
    std::vector<int> pivots( N, 0 );
    // LAPACK information integer
    int info( 0 );
    // Call FORTRAN LAPACK
    LAPACK_DGESV( N, 1, Af.base(), N, &pivots[ 0 ], &( p_B -> operator[] ( 0 ) ), N, info );
    if ( 0 != info )
    {
      std::string problem( " The LAPACK::LU method has detected a failure \n" );
      throw ExceptionExternal( problem, info );
    }
    // compute the determinant sign if asked
    if ( MONITOR_DET )
    {
      DET_SIGN = signature( pivots );
      // product of the diagonals
      for ( std::size_t i = 0; i < N; ++i )
      {
        if ( Af[ i * N + i ] < 0.0 )
        {
          DET_SIGN *= -1;
        }
      }
    }
#endif

  }


  template <>
  void DenseLinearSystem<D_complex>::solve_lapack( )
  {
#ifndef LAPACK
    std::string problem = "The DenseLinearSystem::solve_lapack method has been called\n";
    problem += "but the compiler option -DLAPACK was not provided when\n";
    problem += "the library was built and so LAPACK support is not available.";
    throw ExceptionExternal( problem );
#else

    std::size_t N = p_A -> nrows();
    FortranData Af( *p_A );
    // a vector to store the pivots
    std::vector<int> pivots( N, 0 );
    // LAPACK information integer
    int info( 0 );
    // Call FORTRAN LAPACK
    LAPACK_ZGESV( N, 1, Af.base(), N, &pivots[ 0 ], &( p_B -> operator[]( 0 ).real() ), N, info );
    if ( 0 != info )
    {
      std::string problem( " The LAPACK::LU method has detected a failure \n" );
      throw ExceptionExternal( problem, info );
    }
#endif

  }


  template <typename _Type>
  void DenseLinearSystem<_Type>::solve_native()
  {
    // keep track of row interchanges
    DET_SIGN = 1;
    // first row
    row_iter first_row( p_A -> begin() );
    // last row
    row_iter last_row( p_A -> end() );
    // first RHS elt
    elt_iter first_rhs( p_B -> begin() );
    // last RHS elt
    elt_iter last_rhs( p_B -> end() );
    // iterator start point for the RHS
    elt_iter current_rhs_elt( first_rhs );
    // step through rows
    for ( row_iter current_row = first_row;
          current_row != last_row;
          ++current_row, ++current_rhs_elt )
    {
      // calc an index offset from the first row
      const std::size_t l = distance( first_row, current_row );
      // find max index in column 'l' in the range [l, last_row)
      row_iter max_row = p_A -> max_in_col( l , current_row , last_row );
      // if index of max elt is not the diagonal then swap the rows
      if ( current_row != max_row )
      {
        // flip the determinant sign
        DET_SIGN *= -1;
        // to switch row with max diagonal element to be the current row
        // we could just use
        // std:swap<DenseVector<_Type>>( *current_row, *max_row );
        // but there is no need to swap the eliminated entries
        // to the left of the diagonal
        {
          elt_iter max_row_elt = max_row -> begin() + l;
          for ( elt_iter current_row_elt = current_row -> begin() + l;
                current_row_elt != current_row -> end();
                ++current_row_elt, ++max_row_elt )
          {
            std::swap( *current_row_elt, *max_row_elt );
          }
        }
        // switch the elts in RHS R-vector to be consistent
        std::swap<_Type>( *current_rhs_elt, *( first_rhs + distance( first_row, max_row ) ) );
      }
      // de-reference the element to get diagonal entry value
      const _Type diag_entry = *( current_row -> begin() + l );
      // check it
      if ( std::abs( diag_entry ) < MIN_PIV )
      {
        std::string problem( "A pivot in NDMatrix.GJE is under the minimum tolerance => singular matrix. \n" );
        throw ExceptionRuntime( problem );
      }
      // eliminate all entries below current row in the diagonal's column
      {
        // iterator to the corresponding RHS element
        elt_iter elim_rhs( p_B -> begin() + l + 1 );
        // loop through the rows
        for ( row_iter elim_row = current_row + 1;
              elim_row != last_row;
              ++elim_row, ++elim_rhs )
        {
          // work out the appropriate multiplier for the elimination
          const _Type mult = *( elim_row -> begin() + l ) / diag_entry;
          // elimination of all columns would be *elim_row -= *current_row * mult;
          // but instead optimise out the zero elements to the left of diagonal
          for ( elt_iter col = elim_row -> begin() + l; col != elim_row -> end(); ++col )
          {
            const unsigned col_offset( distance( elim_row -> begin(), col ) );
            *col -= *( current_row -> begin() + col_offset ) * mult;
          }
          // apply the same operation to the RHS vector
          *elim_rhs -= *current_rhs_elt * mult;
        }
      }
    } // close l-loop
    backsub( );
    if ( MONITOR_DET )
    {
      for ( row_iter row = first_row; row != last_row; ++row )
      {
        const unsigned l = distance( first_row, row );
        if ( lt( *( row -> begin() + l ) ) )
        {
          DET_SIGN *= -1;
        }
      }
    }
  }


  template <typename _Type>
  bool DenseLinearSystem<_Type>::lt( _Type value ) const
  {
    std::string problem;
    problem = "You've turned on monitoring of the sign of the determinant for a \n";
    problem += "DenseLinearSystem whose elements are not of type <double>.\n";
    throw ExceptionRuntime( problem );
  }


  template <typename _Type>
  int DenseLinearSystem<_Type>::signature( const std::vector<int> &pivots ) const
  {
    int sign( 1 );
    for ( std::size_t i = 0; i < pivots.size(); ++i )
    {
      if ( pivots[ i ] - 1 != int( i ) )
      {
        sign *= -1;
      }
    }
    return sign;
  }

  template <typename _Type>
  void DenseLinearSystem<_Type >::backsub( ) const
  {
    // last row
    row_riter rend_row( p_A -> rend() );
    // first RHS elt
    elt_riter rbegin_rhs( p_B -> rbegin() );
    // iterator start point for the RHS
    elt_riter current_rhs_elt( rbegin_rhs );
    row_riter current_row( p_A -> rbegin() );
    // step through the rows in reverse order
    for ( ;
          current_row != rend_row;
          ++current_row, ++current_rhs_elt )
    {
      _Type sum = 0.0;
      elt_riter a_elt( current_row -> rbegin() );
      elt_riter b_elt( rbegin_rhs );
      // step through the columns/elts in reverse order
      for ( ;
            b_elt != current_rhs_elt;
            ++a_elt, ++b_elt )
      {
        sum += *a_elt * *b_elt;
      }
      // set the solution using the RHS vector
      *current_rhs_elt = ( *current_rhs_elt - sum ) / *a_elt;
    }
  }


  template class DenseLinearSystem<D_complex>
  ;
  template class DenseLinearSystem<double>
  ;

} // end namespace
