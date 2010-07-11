/// \file Utility.cpp
/// An implementation for a collection of utility functions.

#include <string>
#include <ctime>

#include <Types.h>
#include <Exceptions.h>
#include <Utility.h>
#include <FortranData.h>
#include <FortranBLAS.h>

namespace CppNoddy
{
  namespace Utility
  {

    void time_seed()
    {
      srand( ( unsigned ) std::time( 0 ) );
    }

    void fill_random( DenseVector<double>& V )
    {
      for ( unsigned i = 0; i < V.size(); ++i )
      {
        double x = ( double ) rand() /
                   ( ( double ) RAND_MAX + ( double ) 1 ) ;
        V[ i ] = x;
      }
    }

    void fill_random( DenseVector< D_complex >& V )
    {
      for ( unsigned i = 0; i < V.size(); ++i )
      {
        double x = ( double ) rand() /
                   ( ( double ) RAND_MAX + ( double ) 1 ) ;
        double y = ( double ) rand() /
                   ( ( double ) RAND_MAX + ( double ) 1 ) ;
        V[ i ] = D_complex ( x, y );
      }
    }

    void fill_random( SparseVector<double>& V, const unsigned& num_of_elts )
    {
      do
      {
        double index = ( double ) rand() /
                       ( ( double ) RAND_MAX + ( double ) 1 ) ;
        index *= V.size();
        double x = ( double ) rand() /
                   ( ( double ) RAND_MAX + ( double ) 1 ) ;
        V[ ( unsigned ) index ] = x;
      }
      while ( V.nelts() < num_of_elts );
    }

    void fill_random( SparseVector<D_complex>& V, const unsigned& num_of_elts )
    {
      do
      {
        double index = ( double ) rand() /
                       ( ( double ) RAND_MAX + ( double ) 1 ) ;
        index *= V.size();
        double x = ( double ) rand() /
                   ( ( double ) RAND_MAX + ( double ) 1 ) ;
        double y = ( double ) rand() /
                   ( ( double ) RAND_MAX + ( double ) 1 ) ;
        V[ ( unsigned ) index ] = D_complex ( x, y );
      }
      while ( V.nelts() < num_of_elts );
    }

    void fill_random( DenseMatrix<double>& A )
    {
      DenseVector<double> temp( A.ncols(), 0.0 );
      for ( std::size_t row = 0; row < A.nrows(); ++row )
      {
        fill_random( temp );
        A[ row ] = temp;
      }
    }

    void fill_random( BandedMatrix<double>& A )
    {
      for ( std::size_t row = 0; row < A.nrows(); ++row )
      {
        for ( std::size_t col = std::max( ( int ) ( row - A.noffdiag() ), 0 );
              ( int ) col <= std::min( ( int ) ( row + A.noffdiag() ), ( int ) A.ncols() ); ++col )
        {
          double x = ( double ) rand() / ( ( double ) RAND_MAX + ( double ) 1 ) ;
          A( row, col ) = x;
        }
      }
    }

    DenseVector<double> uniform_node_vector( const double& lower, const double& upper, const std::size_t& N )
    {
      DenseVector<double> V;
      V.reserve( N );
      for ( std::size_t i = 0; i < N; ++i )
      {
        const double delta = ( upper - lower ) / ( N - 1 );
        V.push_back( lower + delta * i );
      }
      return V;
    }

    DenseVector<double> power_node_vector( const double& lower, const double& upper, const std::size_t& N, const double& power )
    {
      DenseVector<double> V( N, 0.0 );
      for ( std::size_t i = 0; i < N; ++i )
      {
        V[ i ] = lower + ( upper - lower ) * std::pow( ( double )i / ( N - 1 ), power );
      }
      return V;
    }

    DenseVector<double> two_uniform_node_vector( const double& lower, const double& mid, const double& upper, const std::size_t& N1, const std::size_t& N2 )
    {
      DenseVector<double> first = uniform_node_vector( lower, mid, N1 );
      DenseVector<double> second = uniform_node_vector( mid, upper, N2 );
      // skip the common elt by starting at 1 not 0
      for ( std::size_t i = 1; i < N2; ++i )
      {
        first.push_back( second[ i ] );
      }
      return first;
    }

    DenseVector<double> mid_weighted_node_vector( const double& lower, const double& upper, const std::size_t& N, const double& weight )
    {
      DenseVector<double> node_vector( N, 0.0 );
      // make a center weighted distribution over -1 to 1
      for ( std::size_t i = 0; i < N; ++i )
      {
        double s( -1.0 + 2.0 * i / ( N - 1 ) );
        node_vector[ i ] = ( weight * std::pow( s, 3 ) + s ) / ( weight + 1 );
      }
      // map the -1 to 1 range to lower to upper
      for ( std::size_t i = 0; i < N; ++i )
      {
        // move the range to 0->2
        node_vector[ i ] += 1.0;
        // move the range tp 0->1
        node_vector[ i ] /= 2.0;
        // move the range to lower -> upper
        node_vector[ i ] *= ( upper - lower );
        node_vector[ i ] += lower;
      }
      return node_vector;
    }

    void vels_from_streamfn_Cartesian( const DenseMatrix<double>& source, const double& dx, const double& dy, DenseMatrix<double>& u, DenseMatrix<double>& v )
    {
      std::size_t Ny( source.nrows() );
      std::size_t Nx( source.ncols() );
      // differentiate the streamfunction to get the velocity field
      {
        // west internal nodes
        std::size_t i( 0 );
        for ( std::size_t j = 1; j < Ny - 1; ++j )
        {
          u( i, j ) = ( source( i, j + 1 ) - source( i, j - 1 ) ) / ( 2 * dy );
          v( i, j ) = -( -source( i + 2, j ) + 4 * source( i + 1, j ) - 3 * source( i, j ) ) / ( 2 * dx );
        }
      }
      {
        // east internal nodes
        std::size_t i( Nx - 1 );
        for ( std::size_t j = 1; j < Ny - 1; ++j )
        {
          u( i, j ) = ( source( i, j + 1 ) - source( i, j - 1 ) ) / ( 2 * dy );
          v( i, j ) = -( source( i - 2, j ) - 4 * source( i - 1, j ) + 3 * source( i, j ) ) / ( 2 * dx );
        }
      }
      {
        // south internal nodes
        std::size_t j( 0 );
        for ( std::size_t i = 1; i < Nx - 1; ++i )
        {
          u( i, j ) = ( -source( i, j + 2 ) + 4 *  source( i, j + 1 ) - 3 * source( i, j ) ) / ( 2 * dy );
          v( i, j ) = -( source( i + 1, j ) - source( i - 1, j ) ) / ( 2 * dx );
        }
      }
      {
        // north internal nodes
        std::size_t j( Ny - 1 );
        for ( std::size_t i = 1; i < Nx - 1; ++i )
        {
          u( i, j ) = ( source( i, j - 2 ) - 4 *  source( i, j - 1 ) + 3 * source( i, j ) ) / ( 2 * dy );
          v( i, j ) = -( source( i + 1, j ) - source( i - 1, j ) ) / ( 2 * dx );
        }
      }
      {
        // corner nodes
        {
          // sw
          std::size_t i( 0 );
          std::size_t j( 0 );
          u( i, j ) = ( -source( i, j + 2 ) + 4 *  source( i, j + 1 ) - 3 * source( i, j ) ) / ( 2 * dy );
          v( i, j ) = -( -source( i + 2, j ) + 4 * source( i + 1, j ) - 3 * source( i, j ) ) / ( 2 * dx );
        }
        {
          // nw
          std::size_t i( 0 );
          std::size_t j( Ny - 1 );
          u( i, j ) = ( source( i, j - 2 ) - 4 *  source( i, j - 1 ) + 3 * source( i, j ) ) / ( 2 * dy );
          v( i, j ) = -( -source( i + 2, j ) + 4 * source( i + 1, j ) - 3 * source( i, j ) ) / ( 2 * dx );
        }
        {
          // ne
          std::size_t i( Nx - 1 );
          std::size_t j( Ny - 1 );
          u( i, j ) = ( source( i, j - 2 ) - 4 *  source( i, j - 1 ) + 3 * source( i, j ) ) / ( 2 * dy );
          v( i, j ) = -( source( i - 2, j ) - 4 * source( i - 1, j ) + 3 * source( i, j ) ) / ( 2 * dx );
        }
        {
          // se
          std::size_t i( Nx - 1 );
          std::size_t j( 0 );
          u( i, j ) = ( -source( i, j + 2 ) + 4 *  source( i, j + 1 ) - 3 * source( i, j ) ) / ( 2 * dy );
          v( i, j ) = -( source( i - 2, j ) - 4 * source( i - 1, j ) + 3 * source( i, j ) ) / ( 2 * dx );
        }
      }
      {
        // interior nodes
        for ( std::size_t i = 1; i < Nx - 1; ++i )
        {
          for ( std::size_t j = 1; j < Ny - 1; ++j )
          {
            u( i, j ) = ( source( i, j + 1 ) - source( i, j - 1 ) ) / ( 2 * dy );
            v( i, j ) = -( source( i + 1, j ) - source( i - 1, j ) ) / ( 2 * dx );
          }
        }
      }
    }

    void vels_from_streamfn_Stokes( const DenseMatrix<double>& source, const double& dr, const double& dz, DenseMatrix<double>& u, DenseMatrix<double>& w )
    {
      std::size_t Nr( source.nrows() );
      std::size_t Nz( source.ncols() );
      // note: this could get a bit iffy if the streamfunction is obtained via a
      // 2nd-order differencing scheme, then we differentiate to get the velocities
      // and have the 1/r factor for small r.
      {
        // east internal nodes
        std::size_t i( Nr - 1 );
        for ( std::size_t j = 1; j < Nz - 1; ++j )
        {
          double r( i * dr );
          u( i, j ) = - ( source( i, j + 1 ) - source( i, j - 1 ) ) / ( 2 * dz ) / r;
          w( i, j ) = ( source( i - 2, j ) - 4 * source( i - 1, j ) + 3 * source( i, j ) ) / ( 2 * dr ) / r;
        }
      }
      {
        // south
        std::size_t j( 0 );
        for ( std::size_t i = 1; i < Nr - 1; ++i )
        {
          double r( i * dr );
          u( i, j ) = - ( -source( i, j + 2 ) + 4 * source( i, j + 1 ) - 3 * source( i, j ) ) / ( 2 * dz ) / r;
          w( i, j ) = ( source( i + 1, j ) - source( i - 1, j ) ) / ( 2 * dr ) / r;
        }
      }
      {
        // north
        std::size_t j( Nz - 1 );
        for ( std::size_t i = 1; i < Nr - 1; ++i )
        {
          double r( i * dr );
          u( i, j ) = - ( source( i, j - 2 ) - 4 * source( i, j - 1 ) + 3 * source( i, j ) ) / ( 2 * dz ) / r;
          w( i, j ) = ( source( i + 1, j ) - source( i - 1, j ) ) / ( 2 * dr ) / r;
        }
      }
      {
        // west a.k.a axis of symmetry
        std::size_t i( 0 );
        for ( std::size_t j = 0; j < Nz; ++j )
        {
          u( i, j ) = 0.0;
          w( i, j ) = 2 * ( source( i + 1, j ) - source( i, j ) ) / ( dr * dr );
        }
        // outer corner nodes
        {
          {
            // ne
            std::size_t i( Nr - 1 );
            std::size_t j( Nz - 1 );
            double r( i * dr );
            u( i, j ) = - ( source( i, j - 2 ) - 4 * source( i, j - 1 ) + 3 * source( i, j ) ) / ( 2 * dz ) / r;
            w( i, j ) = ( source( i - 2, j ) - 4 * source( i - 1, j ) + 3 * source( i, j ) ) / ( 2 * dr ) / r;
          }
          {
            // se
            std::size_t i( Nr - 1 );
            std::size_t j( 0 );
            double r( i * dr );
            u( i, j ) = - ( -source( i, j + 2 ) + 4 * source( i, j + 1 ) - 3 * source( i, j ) ) / ( 2 * dz ) / r;
            w( i, j ) = ( source( i - 2, j ) - 4 * source( i - 1, j ) + 3 * source( i, j ) ) / ( 2 * dr ) / r;
          }
        }
      }
      {
        // interior
        for ( std::size_t i = 1; i < Nr - 1; ++i )
        {
          for ( std::size_t j = 1; j < Nz - 1; ++j )
          {
            double r( i * dr );
            u( i, j ) = - ( source( i, j + 1 ) - source( i, j - 1 ) ) / ( 2 * dz ) / r;
            w( i, j ) = ( source( i + 1, j ) - source( i - 1, j ) ) / ( 2 * dr ) / r;
          }
        }
      }
    }

    DenseVector<double> real( const DenseVector<D_complex>& X )
    {
      DenseVector<double> temp( X.size(), 0.0 );
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        temp[ i ] = X[ i ].real();
      }
      return temp;
    }

    DenseVector<double> imag( const DenseVector<D_complex>& X )
    {
      DenseVector<double> temp( X.size(), 0.0 );
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        temp[ i ] = X[ i ].imag();
      }
      return temp;
    }

    std::string stringify( const int& val )
    {
      std::stringstream temp;
      temp << val;
      return temp.str();
    }

    std::string stringify( const double& val, int p )
    {
      std::stringstream temp;
      temp.precision( p );
      temp << val;
      return temp.str();
    }

    //        // dot product specialised for doubles to use
    //        // the BLAS library if LAPACK specified
    //        double dot<double>( const DenseVector<double>& X, const DenseVector<double>& Y )
    //        {
    //            // check lengths
    //            if ( X.size() != Y.size() )
    //            {
    //                std::string problem( " The LAPACK::dot method has detected a geometry error. \n" );
    //                throw ExceptionGeom( problem, X.size(), Y.size() );
    //            }
    //#ifndef LAPACK
    //            double sum( 0.0 );
    //            inner_product( X.begin(), X.end(), Y.begin(), sum );
    //            return sum;
    //#else
    //            return BLAS_DDOT( X.size(), &X[0], 1, &Y[0], 1 );
    //#endif
    //        }

    DenseMatrix<double> multiply( DenseMatrix<double>& A, DenseMatrix<double>& B )
    {
#ifndef LAPACK
      std::string problem;
      problem = "The Utilities::multiply method has been called\n";
      problem += "but the compiler option -DLAPACK was not provided when\n";
      problem += "the library was built. This non-member function requires BLAS.";
      throw ExceptionExternal( problem );
#else
      // set the matrix geometries
      std::size_t M = A.nrows();
      std::size_t N = B.ncols();
      std::size_t K = A.ncols();
      // No need to transpose first, because we will in fact do B^T * A^T = C^T
      // New the memory for the result
      FortranData Af( A, false );
      FortranData Bf( B, false );
      FortranData Cf( M * N );
#ifdef PARANOID

      if ( K != B.nrows() )
      {
        std::string problem( " The LAPACK::multiply method has detected a failure \n" );
        throw ExceptionGeom( problem, A.nrows(), A.ncols(), B.nrows(), B.ncols() );
      }
#endif
      // call Fortran BLAS
      // call Fortran BLAS
      BLAS_DGEMM ( ( char* ) "N", ( char* ) "N", N, M, K, 1.0, Bf.base(), N, Af.base(), K, 0.0, Cf.base(), N );
      // Return a DenseMatrix<double> from the results -- since Cf is in column_major
      // format, this will actually transpose the Cf data to provide C as required
      return ( Cf.to_dense_matrix( M, N ) );
#endif

    }

  }

} // end namespace
