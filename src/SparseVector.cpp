/// \file SparseVector.cpp
/// Implementation of the SparseVector class -- a sparse, variable size, vector object.
/// The sparse & dense vectors do not share a common base because the
/// sparse class encapsulates an STL map and operator[] assigns entries to
/// the map. Hence get/set methods are used here.

#include <complex>
#include <map>
#include <cmath>
#include <algorithm>

#include <SparseVector.h>
#include <Exceptions.h>
#include <Functors.h>

namespace CppNoddy
{
  template <typename _Type>
  SparseVector<_Type>::SparseVector( const std::size_t& max ) : ZERO( 0.0 ), MAX_SIZE( max )
  {}

  template <typename _Type>
  SparseVector<_Type>::SparseVector( const SparseVector<_Type>& source )
  {
    *this = source;
  }

  template <typename _Type>
  SparseVector<_Type>& SparseVector<_Type>::operator=( const SparseVector<_Type>& source )
  {
    if ( this == &source )
      return * this;
    VEC = source.VEC;
    ZERO = source.ZERO;
    MAX_SIZE = source.MAX_SIZE;
    return *this;
  }

  template <typename _Type>
  SparseVector<_Type>& SparseVector<_Type>::operator*=( const _Type& m )
  {
    for ( iter pos = VEC.begin(); pos != VEC.end(); ++pos )
    {
      pos -> second *= m;
    }
    return *this;
  }

  template <typename _Type>
  SparseVector<_Type>& SparseVector<_Type>::operator-=( const SparseVector<_Type>& X )
  {
#ifdef PARANOID
    if ( X.size() != size() )
    {
      std::string problem;
      problem = " The SparseVector.operator-= method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), X.size() );
    }
#endif
    citer pos_ro = X.VEC.begin();
    iter pos_rw = VEC.begin();
    do
    {
      std::size_t index_rw = pos_rw -> first;
      std::size_t index_ro = pos_ro -> first;
      if ( index_rw == index_ro )
      {
        // element in both vectors
        pos_rw -> second -= pos_ro -> second;
        ++pos_rw;
        ++pos_ro;
      }
      if ( index_rw > index_ro )
      {
        // element is in X but not 'this'
        set
        ( index_ro ) = -( pos_ro -> second );
        ++pos_ro;
      }
      if ( index_rw < index_ro )
      {
        // element is in 'this' but not X
        ++pos_rw;
      }
    }
    while ( pos_ro != X.VEC.end() && pos_rw != VEC.end() );
    if ( pos_ro != X.VEC.end() )
    {
      // need to finish the X data
      do
      {
        set
        ( pos_ro -> first ) = -( pos_ro -> second );
        ++pos_ro;
      }
      while ( pos_ro != X.VEC.end() );
    }
    return *this;
  }

  template <typename _Type>
  SparseVector<_Type>& SparseVector<_Type>::operator+=( const SparseVector<_Type>& X )
  {
#ifdef PARANOID
    if ( X.size() != size() )
    {
      std::string problem;
      problem = " The SparseVector.operator+= method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), X.size() );
    }
#endif
    citer pos_ro = X.VEC.begin();
    iter pos_rw = VEC.begin();
    do
    {
      std::size_t index_rw = pos_rw -> first;
      std::size_t index_ro = pos_ro -> first;
      if ( index_rw == index_ro )
      {
        // element in both vectors
        pos_rw -> second += pos_ro -> second;
        ++pos_rw;
        ++pos_ro;
      }
      if ( index_rw > index_ro )
      {
        // element is in X but not 'this'
        set
        ( index_ro ) = pos_ro -> second;
        ++pos_ro;
      }
      if ( index_rw < index_ro )
      {
        // element is in 'this' but not X
        ++pos_rw;
      }
    }
    while ( pos_ro != X.VEC.end() && pos_rw != VEC.end() );
    if ( pos_ro != X.VEC.end() )
    {
      // need to finish the X data
      do
      {
        set
        ( pos_ro -> first ) = pos_ro -> second;
        ++pos_ro;
      }
      while ( pos_ro != X.VEC.end() );
    }
    return *this;
  }

  template <typename _Type>
  SparseVector<_Type> SparseVector<_Type>::operator+( const SparseVector<_Type>& X ) const
  {
#ifdef PARANOID
    if ( X.size() != size() )
    {
      std::string problem;
      problem = " The SparseVector.operator+ method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), X.size() );
    }
#endif
    SparseVector<_Type> temp( *this );
    temp += X;
    return temp;
  }

  template <typename _Type>
  SparseVector<_Type> SparseVector<_Type>::operator-( const SparseVector<_Type>& X ) const
  {
#ifdef PARANOID
    if ( X.size() != size() )
    {
      std::string problem;
      problem = " The SparseVector.operator- method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), X.size() );
    }
#endif
    SparseVector<_Type> temp( *this );
    temp -= X;
    return temp;
  }

  template <typename _Type>
  SparseVector<_Type> SparseVector<_Type>::operator-() const
  {
    SparseVector<_Type> temp( *this );
    temp *= -1.0;
    return temp;
  }

  template <typename _Type>
  SparseVector<_Type> SparseVector<_Type>::operator*( const _Type& m ) const
  {
    SparseVector<_Type> temp( *this );
    temp *= m;
    return temp;
  }


  template <typename _Type>
  void SparseVector<_Type>::resize( const std::size_t& length )
  {
    MAX_SIZE = length;
  }

  template <typename _Type>
  void SparseVector<_Type>::clear()
  {
    VEC.clear();
  }

  template <typename _Type>
  const _Type SparseVector<_Type>::dot( const SparseVector<_Type>& X ) const
  {
#ifdef PARANOID
    if ( X.size() != size() )
    {
      std::string problem;
      problem = " The SparseVector.dot method is trying to use \n";
      problem += " two vectors of unequal length \n";
      throw ExceptionGeom( problem, size(), X.size() );
    }
#endif
    SparseVector<_Type> temp( X );
    _Type sum( 0.0 );
    for ( iter pos = temp.VEC.begin(); pos != temp.VEC.end(); ++pos )
    {
      std::size_t index = pos -> first;
      /// \todo Using a 'get' (aka find) here is slooow
      /// - obviously this can be improved.
      // the get method returns 0.0 if not stored
      temp[ index ] *= get( index );
    }
    return sum;
  }

  template <typename _Type>
  double SparseVector<_Type>::one_norm() const
  {
    // Accumulate the ABS values of the container entries
    // using a fuction object.
    double sum( 0.0 );
    for ( citer pos = VEC.begin(); pos != VEC.end(); ++pos )
    {
      sum += std::abs( pos -> second );
    }
    return sum;
  }

  template <typename _Type>
  double SparseVector<_Type>::two_norm() const
  {
    // Accumulate the ABS values of the container entries SQUARED
    // using a fuction object.
    double sum( 0.0 );
    for ( citer pos = VEC.begin(); pos != VEC.end(); ++pos )
    {
      sum += std::pow( std::abs( pos -> second ), 2 );
    }
    return sum;
  }

  template <typename _Type>
  double SparseVector<_Type>::inf_norm() const
  {
    double max( 0.0 );
    // Return the maximum (abs) element in the vector
    for ( citer pos = VEC.begin(); pos != VEC.end(); ++pos )
    {
      if ( std::abs( pos -> second ) > max )
      {
        max = std::abs( pos -> second );
      }
    }
    return max;
  }

  template <typename _Type>
  void SparseVector<_Type>::scale( const _Type& scale )
  {
    for ( iter pos = VEC.begin(); pos != VEC.end(); ++pos )
    {
      pos -> second *= scale;
    }
  }

  template <typename _Type>
  void SparseVector<_Type>::add
  ( const SparseVector<_Type>& X )
  {
    *this += X;
  }

  template <typename _Type>
  void SparseVector<_Type>::sub( const SparseVector<_Type>& X )
  {
    *this -= X;
  }

  template <typename _Type>
  std::size_t SparseVector<_Type>::nearest_index( const _Type& value ) const
  {
    citer location( VEC.begin() );
    _Type min_diff( location -> second - value );
    citer start( ++location );
    for ( citer pos = start; pos != VEC.end(); ++pos )
    {
      _Type diff( pos -> second - value );
      if ( std::abs( diff ) < std::abs( min_diff ) )
      {
        min_diff = diff;
        location = pos;
      }
    }
    return location -> first;
  }

  template <typename _Type>
  std::size_t SparseVector<_Type>::maxabs_index() const
  {
    citer location( VEC.begin() );
    double max( std::abs( location -> second ) );
    citer start( ++location );
    for ( citer pos = start; pos != VEC.end(); ++pos )
    {
      if ( std::abs( pos -> second ) > max )
      {
        max = std::abs( pos -> second );
        location = pos;
      }
    }
    return location -> first;
  }

  template <typename _Type>
  std::size_t SparseVector<_Type>::minabs_index() const
  {
    citer location( VEC.begin() );
    double min( std::abs( location -> second ) );
    citer start( ++location );
    for ( citer pos = start; pos != VEC.end(); ++pos )
    {
      if ( std::abs( pos -> second ) < min )
      {
        min = std::abs( pos -> second );
        location = pos;
      }
    }
    return location -> first;
  }

  template <typename _Type>
  void SparseVector<_Type>::swap( const std::size_t& i,
                                  const std::size_t& j )
  {
#ifdef PARANOID
    if ( i >= size() || j >= size() )
    {
      std::string problem;
      problem = " The SparseVector.swap method is trying to access \n";
      problem += " outside the container. \n";
      throw ExceptionRange( problem, size(), std::max( i, j ) );
    }
#endif
    std::swap<_Type>( VEC[ i ], VEC[ j ] );
  }

  template <typename _Type>
  void SparseVector<_Type>::dump() const
  {
    if ( VEC.begin() == VEC.end() )
    {
      std::cout << "[ Empty vector ]\n";
    }
    else
    {
      for ( citer pos = VEC.begin(); pos != VEC.end(); ++pos )
      {
        std::cout << "[ " << pos -> first << " ]" << " = " << pos -> second << " ";
      }
      std::cout << "\n";
    }
  }


  // the templated versions we require are:
  template class SparseVector<double>
  ;
  template class SparseVector<std::complex<double> >
  ;

} // end namespace

